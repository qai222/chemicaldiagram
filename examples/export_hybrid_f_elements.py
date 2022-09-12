import os

os.environ["CSDHOME"] = "/home/qai/local/CCDC/CSD_2022"
from chemicaldiagram import ChemicalDiagram, ChemicalDiagramFilter
import networkx as nx
from chemicaldiagram.converter import from_csd_entry_to_diagram, export_entry
from chemicaldiagram.utils import split_csd_formula, get_identifier_header
from ccdc.io import EntryReader, Entry, csd_version
import re
from tqdm import tqdm
import logging

"""

Task:
---
Search CSD for
```
A list of all f-element + amine containing compounds
- with a very “generous” definition of amine: molecule with a C-N bond, but can have any other types of atoms (O, S, etc.)
```

Before running this:
1. Download and install `CSDS`: send a download request [link](https://www.ccdc.cam.ac.uk/support-and-resources/csdsdownloads/) with your license and activation key.
2. Create a conda env and install CSD API.
    ```
    conda create -n chemicaldiagram
    conda install python==3.7.10  # CSD 2020 is python 3.7
    conda install -y -c <ccdc_conda_offline_channel> Pillow six lxml numpy matplotlib
    conda install -y -c <ccdc_conda_channel> csd-python-api
    # getting werid glibc version error from conda conflict solver if the third command is from <ccdc_conda_channel>
    ```

Possible error on `Ubuntu 20`
```
My PC is:
Description: Ubuntu 20.04.1 LTS
Release: 20.04

Here is the error message shows up after it closes:
libGL error: MESA-LOADER: failed to open i965 (search paths /usr/lib/x86_64-linux-gnu/dri:\$${ORIGIN}/dri:/usr/lib/dri)
libGL error: failed to load driver: i965
libGL error: MESA-LOADER: failed to open i965 (search paths /usr/lib/x86_64-linux-gnu/dri:\$${ORIGIN}/dri:/usr/lib/dri)
libGL error: failed to load driver: i965
libGL error: MESA-LOADER: failed to open swrast (search paths /usr/lib/x86_64-linux-gnu/dri:\$${ORIGIN}/dri:/usr/lib/dri)
libGL error: failed to load driver: swrast
X Error of failed request:  GLXBadContext
  Major opcode of failed request:  152 (GLX)
  Minor opcode of failed request:  6 (X_GLXIsDirect)
  Serial number of failed request:  3215
  Current serial number in output stream:  3214

All backend drivers (.so) are present in the folder ` /usr/lib/x86_64-linux-gnu/dri`. `glxgears` works fine. Even `mercury` works ok.

FWIW, the error message changes if running from an ssh session:
libGL error: MESA-LOADER: failed to open swrast (search paths /usr/lib/x86_64-linux-gnu/dri:\$${ORIGIN}/dri:/usr/lib/dri)
libGL error: failed to load driver: swrast
X Error of failed request:  GLXBadContext
  Major opcode of failed request:  153 (GLX)
  Minor opcode of failed request:  6 (X_GLXIsDirect)
  Serial number of failed request:  3131
  Current serial number in output stream:  3130
```

Solution:
```
rename the file CSD_2020/lib/libstdc++.so.6 to e.g. libstdc++.so.6.bak
```
"""

F6_elements = """La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb"""
F7_elements = """Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No"""
F_elements = set(F6_elements.split() + F7_elements.split())
assert len(F_elements) == 28


def remove_digits(s: str) -> str:
    return re.sub(r'\d+', '', s)


def check_formula(formula: str) -> bool:
    elements = set(re.findall(r"[A-Z]{1}[a-z]{0,1}", formula))

    # the formula should include at least one of the f-elements and C, N
    if not (elements.intersection(F_elements) and {"C", "N"}.issubset(elements)):
        return False

    # the formula should have at least 2 components
    comp_formulae = split_csd_formula(formula)
    if len(comp_formulae) < 2:
        return False

    # at least one component of the formula should include C, N
    comp_element_sets = list(set(remove_digits(f).split()) for f in comp_formulae)
    if not any(cs.issuperset({"C", "N"}) for cs in comp_element_sets):
        return False

    # at least one component of the formula should include F elements
    if not any(cs.intersection(F_elements) for cs in comp_element_sets):
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


def is_f_component(diagram: ChemicalDiagram):
    """
    a diagram is considered as an f component
    1. is connected
    2. contains F element
    3. anion
    """
    assert nx.is_connected(diagram.graph), "chemical checking functions only work for connected graphs"

    if not F_elements.intersection(set(diagram.symbols.values())):
        logging.info("not metal oxide/halide: no M")
        return False

    if diagram.total_charge >= 0:
        logging.info("assumed anion but charge is non-negative")
        return False

    return True


def cdf_function_0(cd: ChemicalDiagram):
    """inclusion: has at least one metal oxide like and at least one amine cation like component"""
    n_f_component = 0
    n_amine_like = 0
    for scd in cd.get_components():
        if is_f_component(scd):
            n_f_component += 1
        # note a component cannot be both amine like and f component
        elif ChemicalDiagramFilter.is_organic_cation_like(scd, is_amine_like=True):
            n_amine_like += 1
    return n_amine_like >= 1 and n_f_component >= 1


def cdf_function_4(cd: ChemicalDiagram):
    """exclusion: less than 2 components"""
    if len(cd.get_components()) < 2:
        return False
    return True


RULES = """line graph inclusion: any C-N-H,"""


def save_json(o, fn):
    with open(fn, 'w') as fp:
        json.dump(o, fp)


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.CRITICAL)
    logging.critical(f"f elements are defined as: {F_elements}")
    CSD_READER = EntryReader('CSD')
    CSD_VERSION = csd_version()

    json_file = "hybrid_f_elements_{}.json".format(CSD_VERSION)
    if os.path.exists(json_file):
        import json

        identifiers = []
        with open(json_file, "r", encoding="utf-8") as f:
            data = json.load(f)
            for d in data:
                identifiers.append(d['diagram']['identifier'])
        with open("hybrid_f_elements_{}.gcd".format(CSD_VERSION), "w", encoding="utf-8") as f:
            for i in identifiers:
                f.write(i + "\n")
    else:
        saved_headers = set()

        in_subgs, ex_subgs = ChemicalDiagramFilter.generate_subgraph_rules(RULES, metal_only=False)

        diagram_filter = ChemicalDiagramFilter(
            in_subgs, ex_subgs,
            [cdf_function_0, cdf_function_4, ]
        )

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
        save_json(save_list, json_file)
