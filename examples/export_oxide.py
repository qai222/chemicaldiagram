import os

os.environ["CSDHOME"] = "/home/csdfordham/CCDC/CSD_2020"

from chemicaldiagram import ChemicalDiagram, ChemicalDiagramFilter
from chemicaldiagram.converter import from_csd_entry_to_diagram, export_entry
from chemicaldiagram.utils import split_csd_formula, MetalAndMetalloids, get_identifier_header
from ccdc.io import EntryReader, Entry
from collections import Counter
import json
from tqdm import tqdm
import re

"""
export first 5 amine-templated metal/metalloid oxides from CSD
"""


def remove_digits(s: str) -> str:
    return re.sub(r'\d+', '', s)


def check_formula(formula: str) -> bool:
    elements = set(re.findall(r"[A-Z]{1}[a-z]{0,1}", formula))
    if not (elements.intersection(MetalAndMetalloids) and {"C", "N", "O"}.issubset(elements)):
        return False
    comp_formulae = split_csd_formula(formula)
    if len(comp_formulae) < 2:
        return False
    comp_element_sets = list(set(remove_digits(f).split()) for f in comp_formulae)
    if not any(cs == {"C", "N", "H"} for cs in comp_element_sets):
        return False
    if not any(cs.intersection(MetalAndMetalloids) for cs in comp_element_sets):
        return False
    return True


def accept_entry(e: Entry, diagram_filter: ChemicalDiagramFilter):
    if not e.has_3d_structure:
        return False
    if not check_formula(e.formula):
        return False
    diagram = from_csd_entry_to_diagram(entry)
    if not diagram_filter.accept(diagram):
        return False
    return True


def cdf_function_0(cd: ChemicalDiagram):
    """inclusion: has at least one metal oxide like and one amine cation like component"""
    n_molike = 0
    n_aminelike = 0
    for scd in cd.get_components():
        if ChemicalDiagramFilter.is_metal_oxide_or_halide(scd, metal_oxide=True, metal_only=False, assume_anion=True):
            n_molike += 1
        if ChemicalDiagramFilter.is_organic_cation_like(scd, is_amine_like=True):
            n_aminelike += 1
    return n_aminelike >= 1 and n_molike >= 1


def cdf_function_1(cd: ChemicalDiagram):
    """inclusion: one of the connected components contains only C N H"""
    sets = []
    for scd in cd.get_components():
        sets.append(set(scd.get_element_list()))
    return {"C", "N", "H"} in sets


def cdf_function_2(cd: ChemicalDiagram):
    """exclusion: any M atom connected to non-oxygen or non-halogen atoms"""
    for k, v in cd.get_env_dict().items():
        if cd.symbols[k] not in MetalAndMetalloids:
            continue
        nb_elements = v["nb_elements"]
        if not set(nb_elements).issubset({"O", "Cl", "F", "Br", "I"}):
            return False
    return True


def cdf_function_3(cd: ChemicalDiagram):
    """exclusion: all M aquo complex"""
    env_dict = cd.get_env_dict()

    def neighboring_hydrogen_count(node):
        return max(
            Counter(env_dict[node]["nb_elements"])["H"],
            Counter(env_dict[node]["nb_elements"])["D"],
        )

    def is_m_aquo_complex(node):
        nb_nodes = env_dict[node]["nb_nodes"]
        nb_nodes = [n for n in nb_nodes if n in cd.local_nodes]  # we only care local oxygen sites
        if any(cd.symbols[n] != "O" for n in nb_nodes):
            return False
        if any(neighboring_hydrogen_count(n) != 2 for n in nb_nodes):
            return False
        return True

    aquo_checks = []
    for k, v in cd.get_env_dict().items():
        if cd.symbols[k] not in MetalAndMetalloids:
            continue
        if is_m_aquo_complex(k):
            aquo_checks.append(True)
        else:
            aquo_checks.append(False)
    if len(aquo_checks) == 0:  # no M at all, excluded
        return False
    if all(aquo_checks):  # all m aquo, excluded
        return False
    else:
        return True


def cdf_function_4(cd: ChemicalDiagram):
    """exclusion: less than 2 components"""
    if len(cd.get_components()) < 2:
        return False
    return True


def cdf_function_5(cd: ChemicalDiagram):
    """exclusion: all M hydroxide"""
    env_dict = cd.get_env_dict()

    def neighboring_hydrogen_count(node):
        return max(
            Counter(env_dict[node]["nb_elements"])["H"],
            Counter(env_dict[node]["nb_elements"])["D"],
        )

    def is_m_hydroxide(node):
        nb_nodes = env_dict[node]["nb_nodes"]
        nb_nodes = [n for n in nb_nodes if n in cd.local_nodes]  # we only care local oxygen sites
        if any(cd.symbols[nb] != "O" for nb in nb_nodes):
            return False
        if not all(neighboring_hydrogen_count(nb) == 1 and len(env_dict[nb]["nb_nodes"]) == 2 for nb in nb_nodes):
            return False
        return True

    hydroxide_checks = []
    for k, v in cd.get_env_dict().items():
        if cd.symbols[k] not in MetalAndMetalloids:
            continue
        if is_m_hydroxide(k):
            hydroxide_checks.append(True)
        else:
            hydroxide_checks.append(False)
    if len(hydroxide_checks) == 0:  # no M at all, excluded
        return False
    if all(hydroxide_checks):  # all m hydroxide, excluded
        return False
    else:
        return True


RULES = """line graph inclusion: any C-N-H,
    line graph exclusion: any O-C-C-C,
    line graph exclusion: any O-C-N,
    line graph exclusion: any O-C-C-N,
    line graph exclusion: any O-C-C-C-N,
    line graph exclusion: any O-C-C-C-C-N,
    line graph exclusion: any P-C-O,
    line graph exclusion: any P-C-C-O,
    line graph exclusion: any P-C-C-C-O,
    line graph exclusion: any O-S-N,
    line graph exclusion: any O-S-C,
    line graph exclusion: any S-C-N,
    line graph exclusion: any P-C-N,
    line graph exclusion: any P-C-C-N,
    line graph exclusion: any P-C-C-C-N,
    line graph exclusion: any P-C-C-C-C-N,
    line graph exclusion: any O-P-C,
    line graph exclusion: any As-C,
    line graph exclusion: any B-C,
    line graph exclusion: any B-N,
    line graph exclusion: any O-N-C,
    line graph exclusion: any O-N-N,
    line graph exclusion: any Si-C,
    line graph exclusion: any M-S-C,
    line graph exclusion: any N-P-O,
    line graph exclusion: any O-P-N,
    line graph exclusion: any O-C-C-O-C-C-O,
    neighbor graph exclusion: any *C has neighbours O;C;C,"""


def save_json(o, fn):
    with open(fn, 'w') as fp:
        json.dump(o, fp)


if __name__ == '__main__':
    # import logging
    # logging.getLogger().setLevel(logging.INFO)

    save_chunk_size = 100

    CSD_READER = EntryReader('CSD')
    saved_headers = set()

    in_subgs, ex_subgs = ChemicalDiagramFilter.generate_subgraph_rules(RULES, metal_only=False)

    diagram_filter = ChemicalDiagramFilter(
        in_subgs, ex_subgs,
        [cdf_function_0, cdf_function_1, cdf_function_2, cdf_function_3, cdf_function_4, cdf_function_5]
    )

    save_list = []
    i_save_list = 0
    for entry in tqdm(CSD_READER):
        identifier = entry.identifier
        header = get_identifier_header(identifier)
        if header in saved_headers:
            continue
        if accept_entry(entry, diagram_filter):
            d = export_entry(entry)
            data = {
                "entry_data": d,
                "diagram": from_csd_entry_to_diagram(entry).as_dict()
            }
            saved_headers.add(header)
            save_list.append(data)
            if len(save_list) == save_chunk_size:
                save_json(save_list, "oxides/{:06}.json".format(i_save_list))
                i_save_list += 1
                save_list = []
