import os

os.environ["CSDHOME"] = "/home/qai/CCDC/CSD_2022"

from chemicaldiagram import ChemicalDiagram, ChemicalDiagramFilter
from chemicaldiagram.converter import from_csd_entry_to_diagram, export_entry
from chemicaldiagram.utils import split_csd_formula, MetalAndMetalloids, get_identifier_header
from ccdc.io import EntryReader, Entry, csd_version
from collections import Counter
import json
from tqdm import tqdm
import re

"""
export entries of organic crystal structures that 1. only one type of molecule 2. has disorder 
"""


def remove_digits(s: str) -> str:
    return re.sub(r'\d+', '', s)


def check_formula(formula: str) -> bool:
    comp_formulae = split_csd_formula(formula)
    if len(comp_formulae) != 1:
        return False
    comp_formula = comp_formulae[0]
    comp_element_set = set(remove_digits(comp_formula).split())
    if not comp_element_set.issuperset({"C", "H"}):
        return False
    if comp_element_set.intersection(MetalAndMetalloids):
        return False
    return True


def accept_entry(e: Entry, diagram_filter: ChemicalDiagramFilter):
    if not e.has_3d_structure:
        return False
    if not check_formula(e.formula):
        return False
    if not e.is_organic:
        return False
    if e.is_organometallic:
        return False
    diagram = from_csd_entry_to_diagram(entry)
    if not diagram_filter.accept(diagram):
        return False
    if not e.has_disorder:
        return False
    return True


def cdf_function_0(cd: ChemicalDiagram):
    """inclusion: no polymetric bonds"""
    return len(cd.graph) == len(cd.local_graph)


def save_json(o, fn):
    with open(fn, 'w') as fp:
        json.dump(o, fp)


if __name__ == '__main__':
    # import logging
    # logging.getLogger().setLevel(logging.INFO)

    CSD_READER = EntryReader('CSD')
    CSD_VERSION = csd_version()
    saved_headers = set()

    diagram_filter = ChemicalDiagramFilter([], [], [cdf_function_0, ])
    size_lim = 100
    save_list = []
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
            if len(save_list) == size_lim:
                break
    save_json(save_list, "disordered_organic/single_organic-{}.json".format(CSD_VERSION))
