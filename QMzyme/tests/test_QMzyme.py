"""
Unit and regression test for the QMzyme package.
"""

# Import package, test suite, and other packages as needed
# Name each function as test_* to be automatically included in test workflow

import sys
import pytest
import QMzyme
import os
from QMzyme.protein_parser import collect_pdb_data, collect_coords

path = os.path.join(os.getcwd(),'etc')
test_start_file='test.pdb'
model = QMzyme.GenerateModel(protein_file=os.path.join(path,test_start_file))
model.catalytic_center(res_name='DNX',res_number=202,chain='A',
                       save_file=True,verbose=False, 
                       output_file='test_catalytic_center.pdb')
model.active_site(distance_cutoff=4,save_file=True,output_file='test_active_site.pdb')
model.truncate(save_file=True, output_file='test_truncate.pdb')

catalytic_center_inp = model.catalytic_center_pdb
active_site_inp = model.active_site_pdb
truncated_active_site_inp = model.truncated_active_site_pdb

catalytic_center_exp = 'test_end_catalytic_center_chainA_DNX_202.pdb'
active_site_exp = 'test_end_catalytic_center_chainA_DNX_202_active_site_distance_cutoff_4.pdb'
truncated_active_site_exp = 'test_end_catalytic_center_chainA_DNX_202_truncated_active_site_distance_cutoff_4.pdb'

def collect_test_data(test_inp,test_exp):
    data_exp = collect_pdb_data(file=os.path.join(path,test_exp),coords=False)
    data_inp = collect_pdb_data(data=test_inp,coords=False)
    #coord_diff = coord_diff_test(data_inp, data_exp)
    return [(data_inp[key],data_exp[key],key) for key in data_exp.keys()]

#def coord_diff_test(inp,exp):
#    a = collect_coords(data=exp)
#    b = collect_coords(data=inp)
#    c = []
#    for i in range(len(a)):
#        if 
#        c.append([a[i][j]-b[i][j] for j in range(3)])
#    return c
        

def test_QMzyme_imported():
    """
    Sample test, will always pass so long as import statement worked.
    """
    assert "QMzyme" in sys.modules

@pytest.mark.parametrize('inp, expected, key',collect_test_data(catalytic_center_inp,catalytic_center_exp))
def test_catalytic_center(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)
    
@pytest.mark.parametrize('inp, expected, key',collect_test_data(active_site_inp,active_site_exp))
def test_active_site(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)

@pytest.mark.parametrize('inp, expected, key',collect_test_data(truncated_active_site_inp,truncated_active_site_exp))
def test_truncate(inp,expected,key):
    assert inp == expected, "a change was detected in {}.".format(key)







