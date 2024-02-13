import QMzyme

model = QMzyme.generate_model(protein_file = 'test_end.pdb')
model.catalytic_center(res_name='DNX',res_number=202,chain='A')
model.active_site(distance_cutoff=4)
model.truncate()

