"""
Unit and regression test for the QMzyme package.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import pytest
import QMzyme
from rdkit import Chem
import os

path = os.path.join(os.getcwd(),'etc')

pdb_format = {'record_type':[0,6],
              'atom_num':[6,11], 
              'atom_name':[12,16],
              'alt_loc':[16],
              'res_name':[17,20],
              'chain_id':[21],
              'res_num':[22,26],
              'insertion_code':[26],
              'x':[30,38],
              'y':[38,46],
              'z':[46,54],
              'occupancy':[54,60],
              'temperature_factor':[60,66],
              'seg_id':[72,76],
              'element_symbol':[76,78],
              'charge':[78,80]}

def pdb_info(info,line):
    return line[pdb_format[info][0]:pdb_format[info][1]]

def grab_data(file):
    try: 
        with open(file, 'r') as f:
            data=f.readlines()
    except:
        raise ValueError("test file {} not found.".format(file))
    
    pdb_data = {'res_name_sequence':[],'atom_sequence':[],'number_of_atoms':0,'atom_coords':[],'res_num_sequence':[]}
    for line in data:
        if pdb_info('record_type',line) in ['ATOM','HETATM']:
            pdb_data['res_name_sequence'].append(pdb_info('res_name',line))
            pdb_data['res_num_sequence'].append(pdb_info('res_num',line))
            pdb_data['number_of_atoms'] += 1
            pdb_data['atom_sequence'].append(pdb_info('atom_name',line))
            pdb_data['atom_coords'].append((pdb_info('x',line),
                                          pdb_info('y',line),
                                          pdb_info('z',line)))

    return pdb_data



def test_QMzyme_imported():
    """
    Sample test, will always pass so long as import statement worked.
    """
    assert "QMzyme" in sys.modules

def test_catalytic_center(main_file='test.pdb', static_test_file='test_catalytic_center_chainA_DNX_202.pdb'):
    """
    tbd
    """
    static_test_file = os.path.join(path,static_test_file)
    static_data = grab_data(static_test_file)
    
    main_file = os.path.join(path,main_file)
    model = QMzyme.generate_model(protein_file=main_file)
    lig = static_test_file.split('_')
    model.catalytic_center(res_name=lig[4],res_number=int(lig[5][:3]),chain=lig[3][-1],save_file=False,output_file='test_cat_center.pdb')
    test_data = grab_data('test_cat_center.pdb')

    print(test_data == static_data)    
    
    
    
    
def test_active_site(test_file='test_catalytic_center_chainA_DNX_202_active_site_distance_cutoff_4.pdb'):
    """
    tbd
    """

def test_truncation(test_file='test_catalytic_center_chainA_DNX_202_truncated_active_site_distance_cutoff_4.pdb'):
    """
    tbd
    """

