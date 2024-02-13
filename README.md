QMzyme
==============================


![](logo.png)


### Jump to Contents
- #### [Installation](#Installation)
- #### [Dependencies](#Dependencies)
- #### [Getting Started](#Getting-Started)
    1. [Generate QMzyme object](#Step-1)
    2. [Define catalytic center](#Step-2)
    3. [Select subsystem](#Step-3)
    4. [Truncate subsystem](#Step-4)
    5. [Resulting model](#Step-5)

## Installation
```
pip install QMzyme
```

## Dependencies
```
- rdkit
- numpy
- aqme
- (optional) MDAnalysis
```

## Getting Started


```python
import QMzyme as qmz
from importlib_resources import files, as_file
import inspect
from IPython.display import HTML, display
```

<a id="Step-1"></a>
## Step 1. Generate your QMzyme object

#### What arguments does QMzyme.GenerateModel take?


```python
print(inspect.signature(qmz.GenerateModel))
```

>     (calculation='QM-only', protein_file=None, pdb_code=None, save_json=True, verbose=True)


#### Initialize model from PDB file:


```python
file = str(files('QMzyme.data').joinpath('1oh0_equ_from_amber_sim.pdb'))
model = qmz.GenerateModel(protein_file=file)

#By default information is printed out after each function call. To turn this off for all calls 
# pass verbose=False to GenerateModel(). Otherwise, you can set verbose True or False in each of the 
# functions individually.
```

>     INITIALIZING... QMZYME OBJECT: 1oh0_equ_from_amber_sim
>     TIMESTAMP: 2024-02-13 09:58:03
    
<a id="Step-2"></a>
## Step 2. Define the catalytic center

#### What arguments does QMzyme.GenerateModel.catalytic_center take?


```python
print(inspect.signature(qmz.GenerateModel.catalytic_center))
```

>     (self, sel='', res_name=None, res_number=None, chain=None, output_file=None, save_file=True, save_json=None, verbose=None)


#### You can specify certain residue identifiers like the residue name and/or residue number. If you use a non-unique specifier (like res_name='ALA'), every Alanine residue will be selected, and a warning will pop up to ensure this was your intention. 


```python
model.catalytic_center(res_number=263, save_file=False)

#I set save_file=False because I don't need the PDB file of the catalytic center,
# but you might want to visualize that to ensure it's what you intended.
```

>     INITIALIZING... CATALYTIC CENTER
>     TIMESTAMP: 2024-02-13 09:58:03.606747
>     DEFINITION:
>     N_ATOMS: 37
>     The following object attributes are now available:
>     self.catalytic_center_definition
>     self.catalytic_center_mol
>     self.catalytic_center_pdb 


#### If you have MDAnalysis installed, you can use the sel argument to pass a string selection command that MDAnalysis will parse:


```python
model.catalytic_center(sel='resid 263', save_json=False, save_file=False)

#I set save_json=False here because an error will pop up if we have 
# the same model object and try to write multiple catalytic centers 
# to its corresponding json file. If you generate a new model object
# in the same directory it will create a new json file and append a
# number identifier to it. 
```

>     INITIALIZING... CATALYTIC CENTER
>     TIMESTAMP: 2024-02-13 09:58:03.822078
>       DEFINITION: 
>       N_ATOMS: 37 
>       The following object attributes are now available:
>    	self.catalytic_center_definition
>     	self.catalytic_center_mol
>     	self.catalytic_center_pdb

<a id="Step-3"></a>
## Step 3. Select the subsystem

#### You must define some cutoff distance in Å. Any residue with at least one atom within this distance of at least one atom of the catalytic center will be included in the subsystem.


```python
model.subsystem(distance_cutoff=5, save_file=False)
```

>     INITIALIZING... SUBSYSTEM SELECTION
>     TIMESTAMP: 2024-02-13 09:58:04.138838
>     CUTOFF: 5
>     N_ATOMS: 427
>     The following object attributes have been generated:
>     self.distance_cutoff
>     self.subsystem_mol
>     self.subsystem_pdb
    

<a id="Step-4"></a>
## Step 4. Truncate the subsystem

#### We will do this following the default truncation scheme of cutting only at backbone atoms of terminal residues (i.e., if the subsystem is comprised of residues 1, 2, 3 and 5, the N termini of residues 1 and 5 will be cut and the C termini of residues 3 and 5 will be cut.


```python
model.truncate(save_file=False)
```

>     INITIALIZING... SUBSYSTEM TRUNCATION
>     TIMESTAMP: 2024-02-13 09:58:04.725904
>     SCHEME: CA_terminal
>     CUTOFF: 5
>     N_ATOMS: 391
>     CHARGE: -1
>     NOTE: charge does NOT include the catalytic center and is based on AMBER amino acid naming conventions.
>     MODEL_COMPONENTS:
>     TYR16,ILE17,VAL20,ASP40,PRO41,TYR57,GLN59,GLY60,LEU61,VAL66,ALA68,MET84,PHE86,VAL88,MET90,LEU99,VAL101,ASH103,MET116,ALA118,TRP120,LEU125,EQU263,WAT372,WAT373,WAT376,WAT378,WAT379,WAT380,WAT385,WAT387,WAT388,WAT389
>     The following object attributes are now available:
>     self.subsystem_charge
>     self.model_atom_count
>     self.truncated_subsystem_mol
>     self.truncated_subsystem_pdb
>     self.constrain_atom_list
>     self.residues
    

<a id="Step-5"></a>
## Step 5. Resulting model

#### As the code was running, a json file was created each step of the way to record information about the creation process and the resulting model. Creation of this json file can be turned off for all steps by setting save_json=False in qmz.GenerateModel(), or at any of the individual steps. Let's take a look at what information is stored in the json:


```python
import json
```


```python
with open(f'{model.protein_prefix}_QMzyme.json') as j:
    model_info = json.load(j)
```


```python
model_info["Starting structure"]
```

>     '/Users/hrk/git/QMzyme/QMzyme/data/1oh0_equ_from_amber_sim.pdb'

```python
model_info["Catalytic center"]
```

>     {'Residue number': 263, 'Number of atoms': 37}




##### You will see there is a section "QMzyme 1" with information on how the subsystem was selected and what residues comprise the truncated subsystem. If you were to run model.subsystem() again, it would create a new entry under "QMyme 2". 


```python
model_info["QMzyme 1"]["Subsystem selection"]
```

>     {'Number of atoms': 427, 'Distance cutoff': 5, 'Output file': 'Not saved'}


```python
model_info["QMzyme 1"]["Truncated subsystem"]
```

>     {'Number of atoms': 391,
>     'Distance cutoff': 5,
>     'Residues': [{'Residue name': 'TYR', 'Residue number': 16, 'Chain': 'A'},
>      {'Residue name': 'ILE', 'Residue number': 17, 'Chain': 'A'},
>      {'Residue name': 'VAL', 'Residue number': 20, 'Chain': 'A'},
>      {'Residue name': 'ASP', 'Residue number': 40, 'Chain': 'A'},
>      {'Residue name': 'PRO', 'Residue number': 41, 'Chain': 'A'},
>      {'Residue name': 'TYR', 'Residue number': 57, 'Chain': 'A'},
>      {'Residue name': 'GLN', 'Residue number': 59, 'Chain': 'A'},
>      {'Residue name': 'GLY', 'Residue number': 60, 'Chain': 'A'},
>      {'Residue name': 'LEU', 'Residue number': 61, 'Chain': 'A'},
>      {'Residue name': 'VAL', 'Residue number': 66, 'Chain': 'A'},
>      {'Residue name': 'ALA', 'Residue number': 68, 'Chain': 'A'},
>      {'Residue name': 'MET', 'Residue number': 84, 'Chain': 'A'},
>      {'Residue name': 'PHE', 'Residue number': 86, 'Chain': 'A'},
>      {'Residue name': 'VAL', 'Residue number': 88, 'Chain': 'A'},
>      {'Residue name': 'MET', 'Residue number': 90, 'Chain': 'A'},
>      {'Residue name': 'LEU', 'Residue number': 99, 'Chain': 'A'},
>      {'Residue name': 'VAL', 'Residue number': 101, 'Chain': 'A'},
>      {'Residue name': 'ASH', 'Residue number': 103, 'Chain': 'A'},
>      {'Residue name': 'MET', 'Residue number': 116, 'Chain': 'A'},
>      {'Residue name': 'ALA', 'Residue number': 118, 'Chain': 'A'},
>      {'Residue name': 'TRP', 'Residue number': 120, 'Chain': 'A'},
>      {'Residue name': 'LEU', 'Residue number': 125, 'Chain': 'A'},
>      {'Residue name': 'EQU', 'Residue number': 263, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 372, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 373, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 376, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 378, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 379, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 380, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 385, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 387, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 388, 'Chain': 'A'},
>      {'Residue name': 'WAT', 'Residue number': 389, 'Chain': 'A'}],
>     'C-alpha atom indices': [2,
>      23,
>      40,
>      54,
>      74,
>      78,
>      97,
>      114,
>      121,
>      138,
>      152,
>      160,
>      175,
>      193,
>      207,
>      222,
>      239,
>      253,
>      264,
>      279,
>      287,
>      309],
>     'Subsystem charge': -1}



#### Copyright
Copyright (c) 2023, Heidi Klem

#### Acknowledgements
Project architecture initiated through the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms).
