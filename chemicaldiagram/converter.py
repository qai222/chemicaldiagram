import datetime
import json
import re
import typing
import pathlib
import inspect

import networkx as nx

from chemicaldiagram import ChemicalDiagram


class ConverterError(Exception): pass


try:
    from ccdc.io import Entry
    from ccdc.io import EntryReader
    from ccdc.search import TextNumericSearch
    from ccdc.entry import Entry, Citation, Molecule, Crystal
    from ccdc.molecule import Bond, Atom
    from ccdc.io import csd_version

    CSD_classes = (Entry, Molecule, Atom, Bond, Crystal, Entry.CrossReference, Citation)
    EntryReader('CSD')
except RuntimeError:
    from os import environ

    if environ.get('CSDHOME') is None:
        error_message = "Environment variable 'CSDHOME' does not exist!"
    else:
        error_message = "'CSDHOME' exists, but failed to import ccdc. Is ccdc installed in the current environment?"
    raise ConverterError(error_message)


class ME(json.JSONEncoder):
    def default(self, o):
        return o.__repr__()


def get_entry_from_identifier(identifier: str) -> Entry:
    csd_reader = EntryReader('CSD')
    return csd_reader.entry(identifier)


def from_csd_entry_to_diagram(entry: Entry) -> ChemicalDiagram:
    identifier = entry.identifier
    e = entry._entry
    cd = e.chemical_diagram()  # csd internal class
    graph = nx.Graph()
    for i in range(cd.natoms()):
        a1 = cd.atom(i)
        side_label = a1.side_label()
        symbol = str(a1.centre_label())
        charge_label = str(a1.top_right_label())

        if len(charge_label) == 0:
            charge = 0
        else:
            try:
                charge = int(re.findall(r"\d+", charge_label)[0])
            except IndexError:
                charge = 1
            if "-" in charge_label:
                charge = - charge
            elif "+" in charge_label:
                pass
            else:
                raise ConverterError("charge_label is strange: {} in {}".format(charge_label, identifier))
        iso_label = str(a1.top_left_label())
        side_element_label = str(side_label[0])
        try:
            side_element_count = int(side_label[1])
        except ValueError:
            if len(side_element_label) > 0:
                side_element_count = 1
            else:
                side_element_count = None
        if side_element_count is None:
            show_label = iso_label + symbol + charge_label + side_element_label
        else:
            show_label = iso_label + symbol + charge_label + side_element_label + str(side_element_count)
        attr = dict(
            symbol=symbol,
            charge_label=charge_label,
            iso_label=iso_label,
            side_element_label=side_element_label,
            side_element_count=side_element_count,
            show_label=show_label,
            charge=charge
        )
        try:
            pos = a1.site().position()
            graph.add_node(i, x=pos.x(), y=pos.y(), **attr)
        except AttributeError:
            graph.add_node(i, x=None, y=None, **attr)
    for i in range(cd.natoms()):
        for j in range(i + 1, cd.natoms()):
            if cd.bond_index(i, j) != -1:
                a1 = cd.atom(i)
                a2 = cd.atom(j)
                b12 = cd.bond(a1, a2)
                bondtype = b12.type().value()
                if b12.polymeric():
                    graph.add_edge(i, j, polymeric=1, bondtype=bondtype)
                else:
                    graph.add_edge(i, j, polymeric=0, bondtype=bondtype)

    # csd by default thinks any bond between O and H is a solid single bond
    # this is problematic if one uses no_all_poly_bond to determine in-box atoms
    # e.g. GITMAX has 3 O where there's at least one non-poly bond
    # solution: if a oxy only poly bond to nH atoms, if it bonds to any H, this bond is changed to poly
    oxynodes_tochange = {}
    symbols = nx.get_node_attributes(graph, "symbol")
    for n in graph:
        if symbols[n] != "O":
            continue
        edge_poly = []
        hnbs = []
        for edgedata in graph.edges(n, data=True):
            n1, n2, eprop = edgedata
            assert n1 == n
            if symbols[n2] in ("H", "D"):
                hnbs.append(n2)
                continue
            edge_poly.append(eprop["polymeric"])
        if all(edge_poly) and len(edge_poly) > 0:
            oxynodes_tochange[n] = hnbs
    for noxy, hnbs in oxynodes_tochange.items():
        for hnb in hnbs:
            graph[noxy][hnb]["polymeric"] = 1
    return ChemicalDiagram(graph, identifier)


def as_dict_crossref(cr: Entry.CrossReference):
    d = {}
    d['type'] = cr.type
    d['text'] = cr.text
    d['identifiers'] = cr.identifiers
    return d


def export_entry(e: Entry) -> dict:
    d = dict()
    d['cif_string'] = e.to_string('cif')
    d['diagram'] = from_csd_entry_to_diagram(e).as_dict()

    attrnames = [a for a in dir(e) if not a.startswith('_')]
    for attrname in attrnames:
        try:
            attr = getattr(e, attrname)
        except RuntimeError:
            continue
        # print(attrname, type(attr))
        if not callable(attr):
            if isinstance(attr, tuple):
                if all(isinstance(item, Entry.CrossReference) for item in attr):
                    d[attrname] = [as_dict_crossref(item) for item in attr]
                    continue
                elif all(isinstance(item, Citation) for item in attr):
                    d[attrname] = [dict(item._asdict()) for item in attr]  # Citation is a named tuple
                    continue
            elif isinstance(attr, datetime.date):
                d[attrname] = attr.isoformat()
            elif type(attr) not in CSD_classes:
                d[attrname] = attr

    j = json.dumps(d, cls=ME)
    d = json.loads(j)
    return d


class ExportedEntry:
    def __init__(
            self,
            diagram: ChemicalDiagram,
            cif_string: str,
            formula: str,
            has_disorder: bool,
            chemical_name: str,
            identifier: str,
    ):
        self.diagram = diagram
        self.cif_string = cif_string
        self.formula = formula
        self.has_disorder = has_disorder
        self.chemical_name = chemical_name
        self.identifier = identifier

    def __repr__(self):
        return "Entry: {}\n\tformula: {}\n\tdiagram: {}".format(self.identifier, self.formula, self.diagram)

    def __hash__(self):
        return hash(self.identifier)

    def __eq__(self, other):
        return self.identifier == other.identifier

    @classmethod
    def from_exported_dict(cls, d):
        diagram = ChemicalDiagram.from_dict(d["diagram"])
        entry_data = d["entry_data"]
        init_dict = {"diagram": diagram}
        signature = inspect.signature(ExportedEntry)
        for p in signature.parameters:
            if p in init_dict:
                continue
            init_dict[p] = entry_data[p]
        return cls(**init_dict)

    def write_cif(self, path: typing.Union[str, pathlib.Path]):
        with open(path, "w") as f:
            f.write(self.cif_string)

    def to_dict(self):
        d = {"diagram": self.diagram.as_dict()}
        signature = inspect.signature(ExportedEntry)
        for p in signature.parameters:
            if p in d:
                continue
            d[p] = getattr(self, p)
        return d

    @classmethod
    def from_dict(cls, d):
        try:
            diagram = ChemicalDiagram.from_dict(d["diagram"])
        except KeyError:
            diagram = None
        init_dict = {"diagram": diagram}

        signature = inspect.signature(ExportedEntry)
        for p in signature.parameters:
            if p in init_dict:
                continue
            init_dict[p] = d[p]

        return cls(**init_dict)
