import logging
from typing import Callable

import networkx as nx

from chemicaldiagram.diagram import ChemicalDiagram
from chemicaldiagram.utils import MetalAndMetalloids, Metals


class ChemicalDiagramFilter:

    def __init__(self,
                 inclusion_subgraphs: [nx.Graph],
                 exclusion_subgraphs: [nx.Graph],
                 cd_filter_functions: [Callable]):

        self.inclusion_subgraphs = inclusion_subgraphs
        self.exclusion_subgraphs = exclusion_subgraphs
        self.cd_filter_functions = cd_filter_functions

    def accept(self, cd: ChemicalDiagram) -> bool:
        for filter_function in self.cd_filter_functions:
            if filter_function(cd):
                logging.info("pass: {}".format(filter_function.__doc__))
            else:
                logging.info("failed: {}".format(filter_function.__doc__))
                return False
        for subg in self.inclusion_subgraphs:
            if cd.contains_subgraph(subg):
                logging.info("pass: {}".format(subg.graph["instruction"]))
            else:
                logging.info("failed: {}".format(subg.graph["instruction"]))
                return False
        for subg in self.exclusion_subgraphs:
            if cd.contains_subgraph(subg):
                logging.info("failed: {}".format(subg.graph["instruction"]))
                return False
            else:
                logging.info("pass: {}".format(subg.graph["instruction"]))
        return True

    @classmethod
    def from_functions_and_rules(cls, filter_functions: [Callable], rules: [str], metal_only: bool = False):
        in_rules, ex_rules = ChemicalDiagramFilter.generate_subgraph_rules(rules, metal_only=metal_only)
        return cls(in_rules, ex_rules, filter_functions)

    @staticmethod
    def gen_line_graph(atoms: [str]):
        g = nx.Graph()
        for i, s in enumerate(atoms):
            g.add_node(i, symbol=s)
        for i in range(len(atoms) - 1):
            g.add_edge(i, i + 1)
        return g

    @staticmethod
    def gen_star_graph(center: str, nbs: [str]):
        """generate a star graph"""
        g = nx.Graph()
        symbols = [center] + nbs
        for i, s in enumerate(symbols):
            g.add_node(i, symbol=s)
        for i in range(1, len(symbols)):
            g.add_edge(0, i)
        return g

    @staticmethod
    def generate_subgraph_rules(rules: str, metal_only: bool = False) -> ([nx.Graph], [nx.Graph]):
        """
        generate inclusion and exclusion rules as (sub)graphs

        :param rules: a string that can be converted to a graph
        :param metal_only: whether `M` means metal, or metal and metalloid
        :return: inclusion graphs, exclusion graphs
        """
        in_subgs = []
        ex_subgs = []

        expand_M_rules = []
        # rules = rules.split("\n")[1:-1]
        rules = rules.split("\n")
        for rule in rules:
            if "M-" in rule:
                if metal_only:
                    for m in Metals:
                        expand_M_rules.append(rule.replace("M-", "{}-".format(m)))
                else:
                    for m in MetalAndMetalloids:
                        expand_M_rules.append(rule.replace("M-", "{}-".format(m)))
        rules = rules + expand_M_rules
        for rule in rules:
            rule = rule.strip().strip(",")
            if "line graph" in rule:
                atoms = [w for w in rule.split() if "-" in w][-1].split("-")
                subg = ChemicalDiagramFilter.gen_line_graph(atoms)
                subg.graph["instruction"] = rule
            elif "neighbor graph" in rule:
                center = [w for w in rule.split() if "*" in w][-1][1:]
                nbs = [w for w in rule.split() if ";" in w][-1].split(";")
                subg = ChemicalDiagramFilter.gen_star_graph(center, nbs)
                subg.graph["instruction"] = rule
            else:
                raise ValueError("rule not understood: {}".format(rule))
            if "inclusion" in rule:
                in_subgs.append(subg)
            elif "exclusion" in rule:
                ex_subgs.append(subg)
            else:
                raise ValueError("rule in/ex not understood: {}".format(rule))
        return in_subgs, ex_subgs

    @staticmethod
    def is_organic_cation_like(diagram: ChemicalDiagram, is_amine_like=False):
        """
        a diagram is considered as an organic cation component if
        1. is connected
        2. contains C H
        3. symbols are allowed elements
        4. total charge positive
        5. node charge non-negative
        6. can be neutralized by removing protons

        in addition, it is considered as an amine cation if
        7. C-N-H is a subgraph
        """
        assert nx.is_connected(diagram.graph), "chemical checking functions only work for connected graphs"

        if not {"C", "H"}.issubset(set(diagram.symbols.values())):
            logging.info("not organic: no coexist carbon and hydrogen")
            return False

        if is_amine_like:
            cnh_graph = ChemicalDiagramFilter.gen_line_graph(["C", "N", "H"])
            if not diagram.contains_subgraph(cnh_graph):
                logging.info("not amine like: no C-N-H substructure")
                return False

        if is_amine_like:
            allowed_symbols = ("C", "N", "H")
        else:
            allowed_symbols = ("C", "N", "H", "F", "Cl", "Br", "I", "S", "O")

        if not set(diagram.symbols.values()).issubset(set(allowed_symbols)):
            logging.info("not organic: elemental composition not a subset of allowed elements")
            return False

        if diagram.total_charge <= 0:
            logging.info("not organic cation as the total charge is non-positive: {}".format(diagram.total_charge))
            return False

        for n, d in diagram.graph.nodes(data=True):
            node_charge = d["charge"]
            if node_charge == 0:
                continue
            if node_charge < 0:
                logging.info("not organic cation as a node has negative charge labelled by CSD: {} {}".format(n, d))
                return False
            else:
                if is_amine_like:
                    if diagram.symbols[n] != "N":
                        logging.info(
                            "amine cation can only have positive node charge on N atoms, but we have: {} {}".format(n,
                                                                                                                    d))
                        return False

                nb_hcount = len([nb for nb in diagram.graph.neighbors(n) if diagram.symbols[nb] == "H"])
                if node_charge > nb_hcount:
                    logging.info("positive atomic charge with fewer protons at: {} {}".format(n, d))
                    return False
        logging.info("pass: this IS an organic cation")
        return True

    @staticmethod
    def is_metal_oxide_or_halide(diagram: ChemicalDiagram, metal_oxide=True, metal_only=False, assume_anion=True):
        """
        a diagram is considered as an metal oxide/halide component if
        1. is connected
        2. contains M

        for oxide:
        - has a MO3 star subgraph

        for halide:
        - has a MX subgraph
        """
        assert nx.is_connected(diagram.graph), "chemical checking functions only work for connected graphs"

        if metal_only:
            m_defined = Metals
        else:
            m_defined = MetalAndMetalloids

        if not set(m_defined).intersection(set(diagram.symbols.values())):
            logging.info("not metal oxide/halide: no M")
            return False

        if assume_anion:
            if diagram.total_charge >= 0:
                logging.info("assumed anion but charge is non-negative")
                return False

        if metal_oxide:
            mo3_subgraphs = []
            for m in m_defined:
                rule = "neighbor graph inclusion: any *{} has neighbours O;O;O,".format(m)
                in_subgs, ex_subgs = ChemicalDiagramFilter.generate_subgraph_rules(rule)
                mo3_subgraphs += in_subgs
            if not any(diagram.contains_subgraph(subg) for subg in mo3_subgraphs):
                logging.info("not metal oxide: no MO3 subgraph")
                return False
        else:
            mx_subgraphs = []
            for m in m_defined:
                for x in ("F", "Cl", "Br", "I"):
                    rule = "line graph exclusion: any {}-{},".format(m, x)
                    in_subgs, ex_subgs = ChemicalDiagramFilter.generate_subgraph_rules(rule)
                    mx_subgraphs += in_subgs
            if not any(diagram.contains_subgraph(subg) for subg in mx_subgraphs):
                logging.info("not metal halide: no MX subgraph")
                return False
        logging.info("pass: this IS a oxide or halide")
        return True
