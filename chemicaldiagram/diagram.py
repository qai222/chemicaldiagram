from __future__ import annotations

import copy
from collections import Counter
from io import StringIO
from operator import eq

import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.isomorphism as iso
import numpy as np

from chemicaldiagram.utils import color_scheme, all_elements

_edges_color_in_drawing = {
    0: ("yellow", "unknown"),
    1: ("black", "single"),
    2: ("red", "double"),
    3: ("blue", "triple"),
    4: ("green", "quadruple"),
    5: ("pink", "aromatic"),
    7: ("brown", "delocalized"),
    9: ("purple", "pi"),
}


class ChemicalDiagram:

    def __init__(self, graph: nx.Graph, identifier: str = None):
        self.graph = graph
        self.identifier = identifier
        self.symbols = nx.get_node_attributes(self.graph, "symbol")

    @property
    def local_graph(self):
        # this is the graph where all nodes are "in the box"
        return self.graph.subgraph(self.get_nodes_not_all_poly()).copy()

    @property
    def local_nodes(self):
        return list(self.local_graph.nodes)

    def __eq__(self, other):
        """edge properties are not considered in determining equality"""
        return hash(self) == hash(other)

    def __len__(self):
        return len(self.local_graph)

    @property
    def graph_hash(self) -> str:
        return nx.weisfeiler_lehman_graph_hash(self.graph, node_attr="symbol", iterations=5, digest_size=64)

    def __hash__(self):
        return hash(self.graph_hash)

    def as_dict(self) -> dict:
        graph = nx.node_link_data(self.graph)
        d = {
            "graph": graph,
            "identifier": self.identifier,
        }
        return d

    @classmethod
    def from_dict(cls, d) -> ChemicalDiagram:
        g = nx.node_link_graph(d["graph"])
        i = d["identifier"]
        return cls(g, i)

    def get_nodes_not_all_poly(self) -> [int]:
        """nodes that connected to at least one non-polymeric bond"""
        ns = []
        for n in self.graph.nodes:
            edge_poly = []
            for edgedata in self.graph.edges(n, data=True):
                n1, n2, eprop = edgedata
                edge_poly.append(eprop["polymeric"])
            if not all(edge_poly):
                ns.append(n)
            elif len(self.graph.edges(n)) == 0:
                ns.append(n)
        return ns

    def get_element_list(self) -> [str]:
        """element list of the local graph"""
        lst = []
        for n in self.local_graph:
            lst.append(self.symbols[n])
        return sorted(lst)

    def get_element_list_all(self) -> [str]:
        lst = []
        for n in self.graph:
            lst.append(self.symbols[n])
        return sorted(lst)

    def get_formula(self):
        local_symbols = []
        for n in self.local_nodes:
            local_symbols.append(self.symbols[n])
        c = dict(Counter(local_symbols))
        s = ""
        for k in sorted(c.keys()):
            s += "{}{}".format(k, c[k])
            s += " "
        return s

    def __repr__(self):
        return "{}: {}; charge: {}.".format(self.__class__.__name__, self.get_formula(), self.total_charge)

    def get_components(self) -> [ChemicalDiagram]:
        components = []
        for c in nx.connected_components(self.graph):
            subgraph = self.graph.subgraph(c).copy()
            components.append(ChemicalDiagram(subgraph, identifier=self.identifier))
        return components

    @property
    def total_charge(self) -> int:
        c = 0
        for n in self.local_graph.nodes(data=True):
            c += n[1]["charge"]
        return c

    def get_env_dict(self, local_nodes_only=True) -> dict:
        """ a lookup table for the elements of a node's neighbors, by default all keys are local nodes"""
        env_dict = {}
        if local_nodes_only:
            nodes = self.local_graph.nodes
        else:
            nodes = self.graph.nodes
        for n in nodes:
            nb_nodes = []
            nb_elements = []
            for nb in self.graph.neighbors(n):
                nb_nodes.append(nb)
                nb_elements.append(self.symbols[nb])
            env_dict[n] = {
                "nb_nodes": sorted(nb_nodes),
                "nb_elements": sorted(nb_elements),
            }
        return env_dict

    def contains_subgraph(self, subgraph) -> bool:
        matcher = iso.GraphMatcher(self.graph, subgraph, node_match=iso.generic_node_match('symbol', None, eq))
        return matcher.subgraph_is_isomorphic()

    def contains_diagram(self, other) -> bool:
        return self.contains_subgraph(other.graph)

    def check_symbols(self):
        if set(self.symbols.keys()).issubset(set(all_elements)):
            return True
        return False

    def draw_svgdata(self, title: str = "", urltxt: str = None, url: str = None):
        f = plt.figure()
        ax = plt.gca()
        ax.set_title(title)
        cdg = self.graph
        posx = nx.get_node_attributes(cdg, 'x')
        posy = nx.get_node_attributes(cdg, 'y')
        cdg_labels = nx.get_node_attributes(cdg, 'show_label')
        cdg_symbols = nx.get_node_attributes(cdg, 'symbol')
        pltgraph = copy.deepcopy(self.graph)
        coords = {}
        subset_symbols = {}
        show_lables = {}
        missingxy = []
        for k in pltgraph.nodes:
            x = posx[k]
            y = posy[k]
            if x is None or y is None:
                missingxy.append(k)
                continue
            coords[k] = (posx[k], posy[k])
            show_lables[k] = cdg_labels[k]
            subset_symbols[k] = cdg_symbols[k]
        for k in missingxy:
            pltgraph.remove_node(k)
        jmolcolors = []
        for n in pltgraph.nodes:
            symb = subset_symbols[n]
            if symb == "D":
                symb = "H"
            jmolcolors.append('#{:02x}{:02x}{:02x}'.format(*color_scheme['Jmol'][symb]))
        nx.draw_networkx_labels(pltgraph, pos=coords, labels=show_lables, ax=ax)
        nx.draw_networkx_nodes(pltgraph, pos=coords, node_color=jmolcolors, ax=ax)

        edge_colors = []
        edge_list = []
        for edge in pltgraph.edges(data=True):
            edge_list.append(edge[:2])
            bt = edge[2]["bondtype"]
            color = _edges_color_in_drawing[bt][0]
            if edge[2]["polymeric"]:
                color = "gray"
            edge_colors.append(color)

        nx.draw_networkx_edges(pltgraph, coords, edge_list, edge_color=edge_colors, ax=ax)

        # nx.draw(pltgraph, with_labels=True, labels=show_lables, pos=coords, ax=ax, node_color=jmolcolors)

        if urltxt and url:
            xycoords_array = [np.array(list(xy)) for xy in coords.values()]
            center = np.mean(xycoords_array, axis=0)
            x, y = center
            plt.text(x, y, urltxt, url=url, bbox=dict(alpha=0.4, url=url, facecolor="red"))

        imgdata = StringIO()
        f.savefig(imgdata, format='svg')
        imgdata.seek(0)  # rewind the data
        svg_dta = imgdata.read()  # this is svg data
        plt.close(f)
        return svg_dta

    def draw_svg(self, title="", urltxt=None, url=None, fn: str = None):
        if fn is None:
            fn = "{}.svg".format(self.identifier)
        data = self.draw_svgdata(title, urltxt, url)
        with open(fn, "w") as f:
            f.write(data)


class BuildingUnitDiagram(ChemicalDiagram):
    _allowed_centers = (
        "Si", "B", "C", "N", "P", "S", "Cl", "As", "Se", "Br", "I",
        "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
        "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba",
        "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
        "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U",
        "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
        "Si", "Te", "Sb", "Ge"
    )

    def __init__(self, g: nx.Graph, identifier=None):
        super().__init__(g, identifier=identifier)

    @property
    def local_graph(self):  # no node is excluded for building units
        return self.graph

    @staticmethod
    def get_pbus_from_diagram(diagram: ChemicalDiagram, allowed_centers=_allowed_centers):
        """
        1. get a list of allowed center nodes
        2. connected components now are chosen as BuildingUnitCores
        3. extend each BuildingUnitCore to its first order neighbours (we can limit such neighbours to be oxygen or hydrogen)

        Resulting subgraph is the BuildingUnit, one for each BuildingUnitCore.
        This is suitable when PBUs share at least one node, such as oxides.
        """
        graph = diagram.graph.copy()
        reduced_graph = diagram.local_graph.copy()
        toberemoved = []
        symbol_dict = diagram.symbols

        for node in reduced_graph.nodes:
            if symbol_dict[node] not in allowed_centers:
                toberemoved.append(node)

        for node in toberemoved:
            reduced_graph.remove_node(node)

        building_units_cores = [reduced_graph.subgraph(c).copy() for c in nx.connected_components(reduced_graph)]
        bus = []
        for buc in building_units_cores:
            building_unit = []
            for buc_node in buc.nodes:
                building_unit.append(buc_node)
                building_unit = building_unit + list(graph.neighbors(buc_node))
            building_unit = list(set(building_unit))
            building_unit_graph = graph.subgraph(building_unit).copy()
            this_bu = BuildingUnitDiagram(building_unit_graph, identifier=diagram.identifier)
            bus.append(this_bu)
        bus = sorted(bus, key=lambda x: len(x), reverse=True)
        return bus
