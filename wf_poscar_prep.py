"""
This script will prepare POSCAR for work-function calculation
it will generate different terminations with the amount of 
vacuum specified
Author : Pedram Tavadze
Email  : petavazohi@mix.wvu.edu
"""

import numpy as np
import pychemia


def remove_atom(st,indices):
    removed_struct = pychemia.core.Structure(symbols=list(np.delete(st.symbols,indices)),cell=st.cell,
                                             positions=list(np.delete(st.positions,indices,axis=0)))
    return removed_struct


def add_layers(structure,nlayer,length,direction):
    st = structure
    alpha = st.lattice.alpha
    beta = st.lattice.beta
    gamma = st.lattice.gamma
    layers = np.zeros(3)
    layers[direction] = nlayer
    a = st.lattice.lengths[0]*(layers[0]+1)
    b = st.lattice.lengths[1]*(layers[1]+1)
    c = st.lattice.lengths[2]*(layers[2]+1)
    newlattice = st.lattice.from_parameters_to_cell(a, b, c, alpha, beta, gamma)
    new_st_added_layers = pychemia.core.Structure(symbols=st.symbols, cell=newlattice.cell, positions=st.positions)
    natom = len(st.positions)
    for ilayer in range(1,nlayer+1):
        for iatom in range(natom) :
            new_st_added_layers.add_atom(st.symbols[iatom],st.positions[iatom]+st.lattice.cell[direction]*ilayer)
    #new_structure = pychemia.core.Structure(symbols=new_structure.symbols, cell=newlattice.cell, positions=new_structure.positions)
    new_st_added_layers.sort_sites()
    n_unique_layers = len(np.unique(st.positions[:,direction].round(decimals=2)))
    structures = []

    for i in range(1,n_unique_layers+1):
        temp_st = new_st_added_layers.copy()
        item = np.unique(temp_st.positions[:,direction])[i*-1]
        layer_index = temp_st.positions[:,direction].round(decimals=2) >= item
        indices = []
        for iatom in range(temp_st.natom):
            if layer_index[iatom] :
                indices.append(iatom)
        temp_st = remove_atom(temp_st,indices=indices)
        a = temp_st.lattice.lengths[0]
        b = temp_st.lattice.lengths[1]
        c = max(temp_st.positions[:,direction][temp_st.positions[:,direction] < item])
        newlattice = st.lattice.from_parameters_to_cell(a, b, c, alpha, beta, gamma)
        new_structure = pychemia.core.Structure(symbols=temp_st.symbols, cell=newlattice.cell, positions=temp_st.positions)
        new_structure = new_structure.add_vacuum(length,direction)
        structures.append(new_structure)
    return structures



st = pychemia.code.vasp.read_poscar("POSCAR")
cs = pychemia.crystal.CrystalSymmetry(st)
st = cs.refine_cell()



structures = add_layers(st,nlayer=2,length=15,direction=2)
i = 0
for new_st in structures :
    pychemia.code.vasp.write_poscar(filepath="POSCAR-with-vac_%i" %i,structure=new_st)
    i += 1


