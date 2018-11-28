"""
This script will prepare POSCAR for work-function calculation
it will generate different terminations with the amount of 
vacuum specified
This script first adds layers then removes atoms from the top according 
to create different terminations
Author : Pedram Tavadze
Email  : petavazohi@mix.wvu.edu
"""

import numpy as np
import pychemia
import argparse

def remove_atom(st,indices):
    structure = pychemia.core.Structure(symbols=list(np.delete(st.symbols,indices)),cell=st.cell,
                                             positions=list(np.delete(st.positions,indices,axis=0)))
    return structure


def add_layers(structure,nlayer,length,direction):
    st = structure
    alpha = st.lattice.alpha
    beta = st.lattice.beta
    gamma = st.lattice.gamma
    n_unique_layers = len(np.unique(st.positions[:,direction]))
    structures = []
    supercell_size=[1,1,1]
    supercell_size[direction] =nlayer
    st = st.supercell(size=(tuple(supercell_size)))
    st.sort_sites()
    structures.append(st.add_vacuum(length,direction))
    for i in range(1,n_unique_layers):
        temp_st = st.copy()
        item = np.unique(temp_st.positions[:,direction])[i*-1]
        layer_index = temp_st.positions[:,direction] >= item
        indices = []
        for iatom in range(temp_st.natom):
            if layer_index[iatom] :
                indices.append(iatom)
        temp_st = remove_atom(temp_st,indices=indices)
        a = temp_st.lattice.lengths[0]
        b = temp_st.lattice.lengths[1]
        c = max(temp_st.positions[:,direction][temp_st.positions[:,direction] < item])
        c = temp_st.lattice.c
        newlattice = st.lattice.from_parameters_to_cell(a, b, c, alpha, beta, gamma)
        new_structure = pychemia.core.Structure(symbols=temp_st.symbols, cell=newlattice.cell, positions=temp_st.positions)
        new_structure = new_structure.add_vacuum(length,direction)
        structures.append(new_structure)
    return structures


if __name__ == "__main__" : 
        
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",dest="input",type=str,help="Input file to add vacuum and layers",default="POSCAR")
    parser.add_argument("--nlayer", dest="nlayer",type=int,help="Number of layers to add",default=3)
    parser.add_argument("--vac_length",dest="length",type=float,help="length of vaccum to add in angstroms",default=15)
    parser.add_argument("--direction",dest="direction",type=int,help="Direction in which the vaccum and layers will be added 0:x ,1:y, 2:z",default=2)
    args = parser.parse_args()
    
    
    infile = args.input
    st = pychemia.code.vasp.read_poscar(infile)
    structures = add_layers(st,nlayer=args.nlayer,length=args.length,direction=args.direction)
    
    i = 0
    for new_st in structures :
        pychemia.code.vasp.write_poscar(filepath="POSCAR-with-vac_%i" %i,structure=new_st,direct=False)
        i += 1
    
    
