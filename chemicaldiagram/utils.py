import json
import os
import re


def load_element_color_scheme():
    module_dir = os.path.dirname(os.path.abspath(__file__))
    json_file = os.path.join(module_dir, "element_color_scheme.json")
    with open(json_file, "r") as f:
        scheme = json.load(f)
    return scheme


color_scheme = load_element_color_scheme()

all_elements = sorted(color_scheme["Jmol"].keys())

Metals = {"Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
          "Ga",
          "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La",
          "Ce",
          "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re",
          "Os",
          "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
          "Cm",
          "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
          "Fl", "Mc", "Lv"}

Metalloids = {"Si", "Te", "Sb", "Ge"}

MetalAndMetalloids = set.union(Metals, Metalloids)


def get_identifier_header(i: str):
    """ return 'ABC' from 'ABC01' """
    return re.findall(r"[A-z]+", i)[0]


def split_csd_formula(formula: str):
    comps = formula.split(",")
    formulae = []
    for comp in comps:
        comp_formula = " ".join(re.findall(r"[A-Z]{1}[a-z]{0,1}\d*", comp))
        if len(comp_formula.strip()) == 0:
            continue
        formulae.append(comp_formula.strip())
    return formulae
