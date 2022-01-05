# Chemical Diagram
Expose `ChemicalDiagram` in `CSD API` as a `networkx.Graph`.

### installation
1. Create a conda env and install `CSD API`.
   ```bash
   conda create -n chemicaldiagram
   conda install python==3.7  # CSD 2020 is python 3.7
   conda install -y -c <ccdc_conda_channel> Pillow six lxml numpy matplotlib
   conda install -y -c <ccdc_conda_channel> csd-python-api
   ```
   Alternatively, you can just use the miniconda environment shipped with CSD-system.
2. Install using `pip`.
   ```bash
   pip install networkx
   pip install chemicaldiagram
   ```
### Usage
```python
from chemicaldiagram.converter import export_entry, from_csd_entry_to_diagram, get_entry_from_identifier
identifier = "GITMAX"
entry = get_entry_from_identifier(identifier)
diagram = from_csd_entry_to_diagram(entry)
print(diagram)
# ChemicalDiagram: C6 Ga1 H21 N2 O10 P2 ; charge: 0.
print(diagram.graph)
# Graph with 47 nodes and 47 edges
```
More examples can be found in the folder of [examples](./examples).
