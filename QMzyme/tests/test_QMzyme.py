"""
Unit and regression test for the QMzyme package.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import pytest
import QMzyme
import os
from QMzyme.utils import collect_pdb_data

path = os.path.join(os.getcwd(),'etc')

def test_QMzyme_imported():
    """
    Sample test, will always pass so long as import statement worked.
    """
    assert "QMzyme" in sys.modules


def catalytic_center_test():
    """
    tbd
    """
    test_start_file='test.pdb'
    test_end_file='test_end_catalytic_center_chainA_DNX_202.pdb'
    
    static_data = collect_pdb_data(os.path.join(path,test_end_file))
    model = QMzyme.generate_model(protein_file=os.path.join(path,test_start_file))
    model.catalytic_center(res_name='DNX',res_number=202,chain='A',
                           save_file=False,verbose=False)
    test_data = collect_pdb_data(data=model.catalytic_center_pdb)
    return [(test_data[key],static_data[key],key) for key in test_data.keys()]

def active_site_test():
    """
    tbd
    """
    test_start_file='test.pdb'
    test_end_file='test_end_catalytic_center_chainA_DNX_202_active_site_distance_cutoff_4.pdb'
    static_data = collect_pdb_data(os.path.join(path,test_end_file))
    
    model = QMzyme.generate_model(protein_file=os.path.join(path,test_start_file))
    model.catalytic_center(res_name='DNX',res_number=202,chain='A',
                           save_file=False,verbose=False)
    model.active_site(distance_cutoff=4,save_file=False)
    
    test_data = collect_pdb_data(data=model.active_site_pdb)
    return [(test_data[key],static_data[key],key) for key in test_data.keys()]

def truncate_test():
    """
    tbd
    """
    test_start_file='test.pdb'
    test_end_file='test_end_catalytic_center_chainA_DNX_202_truncated_active_site_distance_cutoff_4.pdb'
    static_data = collect_pdb_data(os.path.join(path,test_end_file))
    
    model = QMzyme.generate_model(protein_file=os.path.join(path,test_start_file))
    model.catalytic_center(res_name='DNX',res_number=202,chain='A',
                           save_file=False,verbose=False)
    model.active_site(distance_cutoff=4,save_file=False)
    model.truncate(save_file=False)
    
    test_data = collect_pdb_data(data=model.truncated_active_site_pdb)
    return [(test_data[key],static_data[key],key) for key in test_data.keys()]

@pytest.mark.parametrize('inp, expected, key',catalytic_center_test())
def test_catalytic_center(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)
    
@pytest.mark.parametrize('inp, expected, key',active_site_test())
def test_active_site(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)

@pytest.mark.parametrize('inp, expected, key',truncate_test())
def test_truncate(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)




