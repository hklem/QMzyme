QMzyme
==============================


![](logo.png)


### Jump to Contents
- #### [Installation](#Installation)
- #### [Dependencies](#Dependencies)
- #### [Getting Started](#Getting-Started)
    1. [Generate QMzyme object](#Step-1.-Generate-your-QMzyme-object)
    2. [Define catalytic center](#Step-2.-Define-the-catalytic-center)
    3. [Select subsystem](#Step-3.-Select-the-subsystem)
    4. [Truncate subsystem](#Step-4.-Truncate-the-subsystem)
    5. [Resulting model](#Step-5.-Resulting-model)

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
import inspect
```

    /Users/hrk/anaconda3/envs/qmzyme/lib/python3.11/site-packages/MDAnalysis/topology/TPRParser.py:161: DeprecationWarning: 'xdrlib' is deprecated and slated for removal in Python 3.13
      import xdrlib


## Step 1. Generate your QMzyme object

#### What arguments does QMzyme.GenerateModel take?


```python
print(inspect.signature(qmz.GenerateModel))
```

    (calculation='QM-only', protein_file=None, pdb_code=None, save_json=True, verbose=True)


#### Initialize model from PDB file:


```python
model = qmz.GenerateModel(protein_file='QMzyme/tests/init_files/1oh0_equ_xstal.pdb')
# information gets printed out after each function call. To turn this off for all calls 
# pass verbose=False to GenerateModel(). Otherwise, you can turn verbose on or off in each of the 
# functions individually.
```

    INITIALIZING... QMZYME OBJECT: 1oh0_equ_xstal
    TIMESTAMP: 2024-02-12 14:12:27
    


## Step 2. Define the catalytic center

#### What arguments does QMzyme.GenerateModel.catalytic_center take?


```python
print(inspect.signature(qmz.GenerateModel.catalytic_center))
```

    (self, sel='', res_name=None, res_number=None, chain=None, output_file=None, save_file=True, save_json=None, verbose=None)


#### If you have MDAnalysis installed, you use the sel argument to pass a string selection command that MDAnalysis will parse:


```python
model.catalytic_center(sel='resid 263', save_json=False) 
# save_json is false here because an error will pop up if we have 
# the same model object and try to write multiple catalytic centers 
# to its corresponding json file. If you generate a new model object
# in the same directory it will create a new json file and append a
# number identifier to it. 
```

    /Users/hrk/anaconda3/envs/qmzyme/lib/python3.11/site-packages/MDAnalysis/topology/PDBParser.py:331: UserWarning: Element information is missing, elements attribute will not be populated. If needed these can be guessed using MDAnalysis.topology.guessers.
      warnings.warn("Element information is missing, elements attribute "


    INITIALIZING... CATALYTIC CENTER
    TIMESTAMP: 2024-02-12 14:12:27.468466
    DEFINITION: 
    N_ATOMS: 37
    OUTPUT_FILE: 1oh0_equ_xstal_catalytic_center.pdb
    The following object attributes are now available:
    	self.catalytic_center_definition
    	self.catalytic_center_mol
    	self.catalytic_center_pdb


#### Alternatively, you can specify certain components likes the residue name and/or residue number. 


```python
model.catalytic_center(res_number=263)
```

    INITIALIZING... CATALYTIC CENTER
    TIMESTAMP: 2024-02-12 14:12:27.777376
    DEFINITION: 
    N_ATOMS: 37
    OUTPUT_FILE: 1oh0_equ_xstal_catalytic_center.pdb
    The following object attributes are now available:
    	self.catalytic_center_definition
    	self.catalytic_center_mol
    	self.catalytic_center_pdb


## Step 3. Select the subsystem

#### You must define some cutoff distance in Å. Any residue with at least one atom within this distance of at least one atom of the catalytic center will be included in the subsystem.


```python
model.subsystem(distance_cutoff=5)
```

    INITIALIZING... SUBSYSTEM SELECTION
    TIMESTAMP: 2024-02-12 14:12:27.977062
    CUTOFF: 5
    OUTPUT_FILE: 1oh0_equ_xstal_subsystem_distance_cutoff5.pdb
    N_ATOMS: 406
    The following object attributes have been generated:
    	self.distance_cutoff
    	self.subsystem_mol
    	self.subsystem_pdb
    


## Step 4. Truncate the subsystem

#### We will do this following the default truncation scheme of cutting only at backbone atoms of terminal residues (i.e., if the subsystem is comprised of residues 1, 2, 3 and 5, the N termini of residues 1 and 5 will be cut and the C termini of residues 3 and 5 will be cut.


```python
model.truncate()
```

    INITIALIZING... SUBSYSTEM TRUNCATION
    TIMESTAMP: 2024-02-12 14:12:28.518590
    SCHEME: CA_terminal
    CUTOFF: 5
    OUTPUT_FILE: 1oh0_equ_xstal_truncated_subsystem_distance_cutoff5.pdb
    N_ATOMS: 368
    CHARGE: -1
    NOTE: charge does NOT include the catalytic center and is based on AMBER amino acid naming conventions.
     MODEL_COMPONENTS: TYR16,ILE17,VAL20,ASP40,PRO41,TYR57,GLY60,LEU61,VAL66,ALA68,MET84,PHE86,VAL88,MET90,LEU99,VAL101,ASH103,MET105,MET116,ALA118,TRP120,LEU125,EQU263,WAT265,WAT266,WAT267
    The following object attributes are now available:
    	self.subsystem_charge
    	self.model_atom_count
    	self.truncated_subsystem_mol
    	self.truncated_subsystem_pdb
    	self.constrain_atom_list
    	self.residues
    


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




    'QMzyme/tests/init_files/1oh0_equ_xstal.pdb'




```python
model_info["Catalytic center"]
```




    {'Residue number': 263,
     'Number of atoms': 37,
     'Output file': '1oh0_equ_xstal_catalytic_center.pdb'}



##### You will see there is a section "QMzyme 1" with information on how the subsystem was selected and what residues comprise the truncated subsystem. If you were to run model.subsystem() again, it would create a new entry under "QMyme 2". 


```python
model_info["QMzyme 1"]["Subsystem selection"]
```




    {'Number of atoms': 406,
     'Distance cutoff': 5,
     'Output file': '1oh0_equ_xstal_subsystem_distance_cutoff5.pdb'}




```python
model_info["QMzyme 1"]["Truncated subsystem"]
```




    {'Number of atoms': 368,
     'Distance cutoff': 5,
     'Residues': [{'Residue name': 'TYR', 'Residue number': 16, 'Chain': 'A'},
      {'Residue name': 'ILE', 'Residue number': 17, 'Chain': 'A'},
      {'Residue name': 'VAL', 'Residue number': 20, 'Chain': 'A'},
      {'Residue name': 'ASP', 'Residue number': 40, 'Chain': 'A'},
      {'Residue name': 'PRO', 'Residue number': 41, 'Chain': 'A'},
      {'Residue name': 'TYR', 'Residue number': 57, 'Chain': 'A'},
      {'Residue name': 'GLY', 'Residue number': 60, 'Chain': 'A'},
      {'Residue name': 'LEU', 'Residue number': 61, 'Chain': 'A'},
      {'Residue name': 'VAL', 'Residue number': 66, 'Chain': 'A'},
      {'Residue name': 'ALA', 'Residue number': 68, 'Chain': 'A'},
      {'Residue name': 'MET', 'Residue number': 84, 'Chain': 'A'},
      {'Residue name': 'PHE', 'Residue number': 86, 'Chain': 'A'},
      {'Residue name': 'VAL', 'Residue number': 88, 'Chain': 'A'},
      {'Residue name': 'MET', 'Residue number': 90, 'Chain': 'A'},
      {'Residue name': 'LEU', 'Residue number': 99, 'Chain': 'A'},
      {'Residue name': 'VAL', 'Residue number': 101, 'Chain': 'A'},
      {'Residue name': 'ASH', 'Residue number': 103, 'Chain': 'A'},
      {'Residue name': 'MET', 'Residue number': 105, 'Chain': 'A'},
      {'Residue name': 'MET', 'Residue number': 116, 'Chain': 'A'},
      {'Residue name': 'ALA', 'Residue number': 118, 'Chain': 'A'},
      {'Residue name': 'TRP', 'Residue number': 120, 'Chain': 'A'},
      {'Residue name': 'LEU', 'Residue number': 125, 'Chain': 'A'},
      {'Residue name': 'EQU', 'Residue number': 263, 'Chain': 'A'},
      {'Residue name': 'WAT', 'Residue number': 265, 'Chain': 'A'},
      {'Residue name': 'WAT', 'Residue number': 266, 'Chain': 'A'},
      {'Residue name': 'WAT', 'Residue number': 267, 'Chain': 'A'}],
     'C-alpha atom indices': [2,
      23,
      40,
      54,
      74,
      78,
      97,
      104,
      121,
      135,
      143,
      158,
      176,
      190,
      205,
      222,
      236,
      247,
      262,
      277,
      285,
      307],
     'Subsystem charge': -1,
     'Output file': '1oh0_equ_xstal_truncated_subsystem_distance_cutoff5.pdb'}



#### Copyright
Copyright (c) 2023, Heidi Klem

#### Acknowledgements
Project architecture initiated through the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms).
