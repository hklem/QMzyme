


##### PDB SPECIFIC #####
col_format = {'record_type':[0,6],
              'atom_number':[6,11], 
              'atom_name':[12,16],
              'alt_loc':[16],
              'res_name':[17,20],
              'chain_id':[21],
              'res_number':[22,26],
              'insertion_code':[26],
              'x':[30,38],
              'y':[38,46],
              'z':[46,54],
              'occupancy':[54,60],
              'temperature_factor':[60,66],
              'seg_id':[72,76],
              'element_symbol':[76,78],
              'charge':[78,80]}

pdb_format = {
    'record_type': lambda line: line[0:6].split()[0],
    'atom_number': lambda line: int(line[6:11].split()[0]), 
    'atom_name': lambda line: line[12:16].split()[0],
    'alt_loc': lambda line: line[16].split()[0],
    'res_name': lambda line: line[17:20].split()[0],
    'chain_id': lambda line: line[21].split()[0],
    'res_number': lambda line: int(line[22:26].split()[0]),
    'insertion_code': lambda line: line[26].split()[0],
    'atom_coords': lambda line: tuple([float(coord) for coord in line[30:54].split()]),
    'occupancy': lambda line: line[54:60].split()[0],
    'temperature_factor': lambda line: line[60:66].split([0]),
    'seg_id': lambda line: line[72:76].split()[0],
    'element_symbol': lambda line: line[76:78].split()[0],
    'charge': lambda line: int(line[78:80].split()[0]),
    }
col_info_map = {
    'HEADER': ['CLASSIFICATION', 'DEPOSIT DATE', 'PDB ID'],
    'SOURCE': ['ORGANISM_SCIENTIFIC'],
    'EXPDTA': ['EXPDTA'],
    'AUTHOR': ['AUTHOR'],
    'JRNL': ['DOI'],
    'REMARK': ['RESOLUTION.','TEMPERATURE',' PH ','MUTATION'],
    'HET ': ['HET'],
    'HETNAM': ['HETNAM']
    }

def pdb_deposit_info(file):
    deposit_info = {}
    with open(file) as f:
        data = f.readlines()
    deposit_info['HET'] = []
    deposit_info['HETNAM'] = []
    deposit_info['MUTATION'] = []
    for line in data:
        for key in col_info_map.keys():
            if line.startswith(key):
                for i,info in enumerate(col_info_map[key]):
                    if key == 'HEADER':
                        deposit_info[info] = line.split(key)[1].split()[i]
                    elif key in ['HET ', 'HETNAM']:
                        deposit_info[info].append(line.split(info)[1].split())
                    elif info in line:
                        if info == 'MUTATION':
                            deposit_info[info].append(line.split(info)[-1].split())
                        else:
                            try: deposit_info[info]
                            except: deposit_info[info] = line.split(info)[1]
                            
    for item in col_info_map.items():
        for key in item[1]:
            try: 
                deposit_info[key]
                if deposit_info[key] == []:
                    deposit_info[key] = None
            except: deposit_info[key] = None
            
    deposit_info['TEMPERATURE (K)'] = float(deposit_info['TEMPERATURE'].split(':')[1])
    deposit_info['PH'] = deposit_info[' PH '].split(':')[-1].split()[0]
    deposit_info['ORGANISM_SCIENTIFIC'] = deposit_info['ORGANISM_SCIENTIFIC'].split(':')[1]
    deposit_info['ORGANISM_SCIENTIFIC'] = deposit_info['ORGANISM_SCIENTIFIC'].split(';')[0]
    deposit_info['RESOLUTION (Å)'] = float(deposit_info['RESOLUTION.'].split()[0])
    del deposit_info['RESOLUTION.']
    del deposit_info['TEMPERATURE']
    del deposit_info[' PH ']
        
    #deposit_info['RESOLUTION'] = 
    return deposit_info

        #try:
        #    if line.startswith('HEADER'):
        #        deposit_info['PDB ID'] = line.split()[-1]
        #        deposit_info['CLASSIFICATION'] = line.split()[1]
        #        deposit_info['DEPOSIT DATE'] = line.split()[2] 
        #    elif 'ORGANISM' in line.split()[2]:
        #        deposit_info[line.split()[2].split(':')[0]] = line.split()[3:]
        #    elif line.startswith('AUTHOR'):
        #        deposit_info['AUTHOR'] = line.split()[1:]
        #    elif line.startswith('EXPDTA'):
        #        deposit_info['EXPERIMENT'] = line.split()[1:]
        #    elif 'RESOLUTION.' in line and line.startswith('REMARK'):
        #        try: 
        #            deposit_info['RESOLUTION']
        #        except:
        #            deposit_info['RESOLUTION'] = float(line.split('RESOLUTION.')[1].split()[0])
        #    elif 'DOI' in line:
        #        deposit_info['DOI'] = line.split('DOI')[1]
        #    elif line.split()[0] =='HET':
        #        het.append(line.split('HET')[-1])
        #    elif line.split()[0] == 'HETNAME':
        #        het_name.append(line.split('HETNAM')[-1])
        #except:
        #    pass
    #deposit_info['NON PROTEIN COMPONENTS'] = het
    #deposit_info['NON PROTEIN NAMES'] = het_name
    #return deposit_info

def pdb_info(line, info):
    return pdb_format[info](line)

def find_file_or_data(file, data):
    if file is not None:
        try: 
            with open(file, 'r') as f:
                data=f.readlines()
            return data
        except:
            raise FileNotFoundError("file {} not found.".format(file))
    elif data is None:
        raise Exception("Must define either file or lines from a pdb file.")
    else:
        return data

def collect_pdb_data(file=None, 
                     data=None, 
                     info=list(pdb_format.keys()), 
                     verbose=False):
    data = find_file_or_data(file,data)
    pdb_data = {}
    for key in info:
        pdb_data[key] = []
    for line in data:
        if pdb_format['record_type'](line) in ['ATOM','HETATM']:
            for key in info:
                try: pdb_data[key].append(pdb_format[key](line))
                except: pdb_data[key].append(None)
    for key in info:
        try:
            if list(set(pdb_data[key])) == [None]:
                pdb_data[key] = None
                if verbose is True:
                    print("WARNING: No {} information found in pdb data.".format(key))
        except: pass
    return pdb_data

