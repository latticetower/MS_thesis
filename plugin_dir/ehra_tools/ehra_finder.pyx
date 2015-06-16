#!/usr/bin/env python
from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

from benchmarkers import *
from utils import *
from data_readers import read_pdb_info, read_dssp_info, process_chain

import pyximport; pyximport.install()

from ehra_interface import *
from ehra_pockets import *
from ehra_loops import *

#input: array of triangles
#output: aminoacids
def to_aa(triangles, surface):
    if len(triangles) == 0:
        return triangles
    g = np.vectorize(lambda x: surface[0][x].get_parent())
    return np.unique(g(triangles))

@timed
def find_regions(pdb_filename, chain1, chain2, cutoff = 5.0, sas_radius = 1.4, masked_aminoacids = []):
    pair_of_chains = read_pdb_info(pdb_filename, chain1, chain2, masked_aminoacids)
    if pair_of_chains == None:
        return None
    chains_ss_info = read_dssp_info(pdb_filename, chain1, chain2)
    if chains_ss_info == None:
        return None
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    triangles = extend_interface_1(triangles, surface1)
    ## custom function: returns True if can travel from one triangle to another
    def check_edge(p1, p2):
        x1 = surface1[0][p1]
        x2 = surface1[0][p2]
        dist = (x1 - x2) - (get_radius(x1.element) + get_radius(x2.element) + sas_radius*2)
        return dist > 0
    def check_simplex(triangle):
        from alpha_shapes import check_triangle
        points = np.vectorize(lambda x: surface1[0][x])(triangle)
        return check_triangle(*points, r=sas_radius)
    def distance_func(x1, x2):
        return x1 - x2
    G = DTGraph(surface1)
    pockets = G.find_pockets(triangles, check_edge, distance_func, check_simplex)
    #print "pockets {0}".format(pockets)
    nodes = diff(pockets, triangles)
    #print "nodes {0}".format(nodes)
    #print(nodes)
    aa_with_coils = extend_to_coils(pockets, chains_ss_info[chain1], surface1)
    chain_length = pair_of_chains[2]
    nn1 = np.union1d(
        to_aa(aa_with_coils, surface1),
        to_aa(pockets, surface1))
    nn2 = to_aa(triangles, surface1)
    #print aa_with_coils
    #    if (len(nodes) > 0):
    if len(nn1) == 0:
        if len(nn2) == 0:
            return (nn1, chain_length)
        return (nn2, chain_length)
    else:
        if len(nn2) == 0:
            return (nn1, chain_length)
    return (np.unique(np.union1d(
        nn1,
        nn2
        )), #this retuns all aminoacids forming pockets except ones shown previously
        chain_length
    )
    #C = CHGraph(surface1)
