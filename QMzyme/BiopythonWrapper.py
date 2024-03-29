import os
import numpy as np
from typing import List, Optional
from QMzyme.Biopython.StructureBuilder import StructureBuilder
from QMzyme.Biopython.Structure import Structure
from QMzyme.Biopython import NeighborSearch
from QMzyme.Biopython.PDBParser import PDBParser
from QMzyme.Biopython.MMCIFParser import MMCIFParser
from QMzyme.Biopython.PDBIO import PDBIO
from QMzyme.Biopython.PDBIO import StructureIO
from QMzyme.utils import filename_format

model_dispatch = {'id': lambda model: model.id,
                  'full_id': lambda model: model.full_id,
                  'serial_number': lambda model: model.serial_number,
                  'chains_list': lambda model: model.child_list,
                  'chains_dict': lambda model: model.child_dict,
                  'xtra': lambda model: model.extra,
                  'structure': lambda model: model.parent}

chain_dispatch = {'id': lambda chain: chain.id,
                  'full_id': lambda chain: chain.full_id,
                  'residues_list': lambda chain: chain.child_list,
                  'residues_dict': lambda chain: chain.child_dict,
                  'xtra': lambda chain: chain.extra,
                  'model': lambda chain: chain.parent,
                  'structure': lambda chain: model_dispatch['structure'](chain.parent)}

res_dispatch = {'resname': lambda res: res.resname,
                'resnumber': lambda res: res.id[1],
                'chain': lambda res: res.parent.id,
                'segid': lambda res: res.segid,
                'id': lambda res: res.id,
                'full_id': lambda res: res.full_id,
                'atoms_list': lambda res: res.child_list,
                'atoms_dict': lambda res: res.child_dict,
                'xtra': lambda res: res.extra,
                'model': lambda res: chain_dispatch['model'](res.parent),
                'structure': lambda res: chain_dispatch['structure'](res.parent)}

atom_dispatch = {'name': lambda atom: atom.name,
                 'fullname': lambda atom: atom.fullname,
                 'coord': lambda atom: atom.coord,
                 'b_factor': lambda atom: atom.bfactor,
                 'bfactor': lambda atom: atom.bfactor,
                 'occupancy': lambda atom: atom.occupancy,
                 'altloc': lambda atom: atom.altloc,
                 'full_id': lambda atom: atom.full_id,
                 'id': lambda atom: atom.id,
                 'serial_number': lambda atom: atom.serial_number,
                 'element': lambda atom: atom.element,
                 'mass': lambda atom: atom.mass,
                 'pqr_charge': lambda atom: atom.pqr_charge,
                 'radius': lambda atom: atom.radius,
                 'xtra': lambda atom: atom.extra,
                 'residue': lambda atom: atom.parent,
                 'resname': lambda atom: res_dispatch['resname'](atom.parent),
                 'resnumber': lambda atom: res_dispatch['resnumber'](atom.parent),
                 'chain': lambda atom: res_dispatch['chain'](atom.parent),
                 'segid': lambda atom: res_dispatch['segid'](atom.parent),
                 'model': lambda atom: res_dispatch['model'](atom.parent),
                 'structure': lambda atom: res_dispatch['structure'](atom.parent)}

class BiopythonWrapper():

    def write_pdb(pdb_object, filename: Optional[str]=""):
        io = PDBIO()
        io.set_structure(pdb_object)
        if filename == "":
            for i in pdb_object.full_id:
                filename += f"{i}_"
            filename = filename[:-1]
        filename = filename_format(filename,'.pdb')
        io.save(f"{filename}")
        print(f"File {filename} created.")

    def order_residues(residues_list):
        new_list = []
        chains = list(set([res.parent.id for res in residues_list]))
        for chain in chains:
            res_list = [res for res in residues_list if res.parent.id == chain]
            max = np.max([res.id[1] for res in res_list])
            min = np.min([res.id[1] for res in res_list])
            for i in range(min, max+1):
                for res in res_list:
                    if res.id[1] == i:
                        new_list.append(res)
        return new_list

    def count_atoms(entity):
        return len(list(entity.get_atoms()))

    def count_residues(entity):
        return len(list(entity.get_residues()))

    def count_chains(entity):
        return len(list(entity.get_chains()))
    
    def load_structure(file, id):
        if file.endswith('pdb'):
            p = PDBParser(QUIET=True)
        elif file.endswith('cif'):
            p = MMCIFParser(QUIET=True)
        else:
            print(f"File format for {file} not recognized. Please use filename ending in '.pdb' or '.cif'.")
        if id is None:
            id = os.path.basename(file).split('.')[0]
        structure = p.get_structure(f'{id}', f'{file}')
        return structure
    
    def get_neighbors(entity, coords: [list], cutoff):
        #return will be one list of all atoms that were within cutoff
        neighbors = []
        neighbor_search = NeighborSearch(list(entity.get_atoms()))
        if len(np.shape(coords)) == 1:
            coords = [coords]
        for coord in coords:
            neighbors += neighbor_search.search(coord,cutoff)
        return list(set(neighbors)) 
    
    def init_model(structure, model_residues, xtra = {}):
        """
        xtra is a dictional where you can add whatever you want.
        """ 
        s = StructureBuilder()
        s.init_structure(structure.id)
        s.init_model(model_id = len(list(structure.get_models())))
        for c in structure.child_list[0].get_chains():
            #for res in self.models[-1]:
            for res in model_residues:
                if res.parent.id != c.id:
                    continue
                s.init_chain(chain_id = c.id)        
                s.init_seg(segid = ' ')
                resname = res.resname
                field = res.id[0]
                resseq = res.id[1]
                icode = res.id[2]
                s.init_residue(resname, field, resseq, icode)
                for atom in res.get_atoms():
                    name = atom.name
                    coord = atom.coord
                    b_factor = atom.bfactor
                    occupancy = atom.occupancy
                    altloc = atom.altloc
                    fullname = atom.fullname
                    serial_number = atom.serial_number
                    element = atom.element
                    pqr_charge = atom.pqr_charge
                    radius = atom.radius
                    is_pqr = False
                    if pqr_charge is not None:
                        is_pqr = True
                    s.init_atom(name, coord, b_factor, occupancy, altloc, fullname, serial_number, element, pqr_charge, radius, is_pqr)
        s.model.xtra = xtra

        return s.model








