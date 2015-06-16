#!/usr/bin/env python

from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

from benchmarkers import *
from utils import *

def get_chain_atoms_without_water_and_masked_aa(structure, chain, masked_aminoacids = []):
    return [
        atom for atom in structure[0][chain].get_atoms()
            if atom.get_parent().get_id()[0] != 'W' and
                not atom.get_parent().get_id()[1] in masked_aminoacids
        ]

@timed
def read_pdb_info(filename, chain1= 'L', chain2 = 'H', masked_aminoacids = []):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    if chain2 == None:
        chain_names = map(lambda x: x.id, structure[0])
        if len(chain_names) != 2:
            print "there are more than 2 chains, can't decide"
            return None
        if chain1 != None:
            chain2 = chain_names[1] if chain1 == chain_names[0] else chain_names[0]
    chain1_data = get_chain_atoms_without_water_and_masked_aa(structure, chain1, masked_aminoacids)
    chain2_data = get_chain_atoms_without_water_and_masked_aa(structure, chain2)
    chain1_residues = {atom.get_parent().get_id()[1] for atom in chain1_data}
    #print chain1_residues
    return (
        chain1_data,
        chain2_data,
        len(chain1_residues)
    )

@timed
def read_pdb_chains_info(filename):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    return [ chain.get_id() for chain in structure[0] ]

def group_elements(data):
    from itertools import groupby
    coils = ['-', 'S', 'T', 'B', 'C', 'G']
    last_coil = None
    start_coil_pos = -1
    protein_coils = [] #contains tuples - (start position, end upper bound position, coil type)
    last_key = -1
    #print data
    l = [el for el in data if el[1] in coils]
    protein_coils = []
    for (k, v) in l:
        if last_key != k - 1:
            if last_key > 0:
                protein_coils.append((start_coil_pos, last_key))
            start_coil_pos = k
        last_key = k
    if last_key > 0:
        protein_coils.append((start_coil_pos, last_key))
    return protein_coils

@timed
def read_dssp_info(filename,
                    chain1= 'L', chain2 = 'H',
                    key_transf1 = lambda x: (x.get_id()[1]),
                    key_transf2 = lambda x: x.get_parent().get_id()):
    chains_data = {}
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    model = structure[0]
    dssp = DSSP(model, filename, dssp='mkdssp')
    if chain2 == None:
        chain_names = map(lambda x: x.id, structure[0])
        if len(chain_names) != 2:
            print "there are more than 2 chains, can't decide"
            return None
        if chain1 != None:
            chain2 = chain_names[1] if chain1 == chain_names[0] else chain_names[0]
    for chain in (chain1, chain2):
        data = sorted([
                (key_transf1(key), value)
                for (key, value, v, e, r,u) in dssp
                if key_transf2(key) == chain
            ])
        chains_data[chain] = group_elements(data)
    return chains_data

# second arg - to test on very-very small subset
# should return triangulation, with some additional objects - neighbours of edge and neighbours of vertex
@timed
def process_chain(points, points_no = 0):
    print("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay, KDTree
    p0 = points[: (points_no if points_no > 0 else len(points))]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    tri = Delaunay(p)
    tree = KDTree(p)
    return (p0, tri, tree, p)
