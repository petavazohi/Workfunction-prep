# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 17:19:00 2018

@author: Pedram Tavadze
"""
import os
import pychemia

ls = os.listdir('.')
direction = 2 # z 
for item in ls : 
    if "POSCAR" in item : 
        st = pychemia.code.vasp.read_poscar(item)
        print (item,st.lattice.c-max(st.positions[:,direction])) 