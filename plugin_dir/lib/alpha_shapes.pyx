#!/usr/bin/env python

from Bio.PDB import *
import numpy as np
import timeit, logging, os, math

from benchmarkers import *
from utils import get_radius

# this scripts gets alpha shape for given triangulated set of atoms from pdb
# as a union of balls
# uses qhull.org python interface from scipy

def volume_2(l1_2, l2_2, l3_2, l4_2, l5_2, l6_2):
    return (l1_2 * l5_2 * (l2_2 + l3_2 + l4_2 + l6_2 - l1_2 - l5_2) +
     l2_2 * l6_2 * (l1_2 + l3_2 + l4_2 + l5_2 - l2_2 - l6_2) +
     l3_2 * l4_2 * (l1_2 + l2_2 + l5_2 + l6_2 - l3_2 - l4_2) -
     l1_2 * l2_2 * l4_2 -
     l2_2 * l3_2 * l5_2 -
     l1_2 * l3_2 * l6_2 -
     l4_2 * l5_2 * l6_2
    ) / 144.0
def area_2(a_2, b_2, c_2):
    result = (
     2 * a_2 * b_2 + 2 * a_2 * c_2 + 2 * b_2 * c_2 - a_2 * a_2 - b_2 * b_2 - c_2 * c_2
    ) / 16.0
    if result >= 0.0:
        return math.sqrt(result)
    return -1.0

# calls get_radius(atom_name)
# gets 4 atoms
# returns True if difference between their van der waals spheres and tetrahedra built on their atom centers is zero,
# False - otherwise.
def check_names(a1, a2, a3, a4, coords_method, name_method):
    r1 = get_radius(a1.element) ** 2
    r2 = get_radius(a2.element) ** 2
    r3 = get_radius(a3.element) ** 2
    r4 = get_radius(a4.element) ** 2
    l1 = (a1 - a2) ** 2
    l2 = (a1 - a3) ** 2
    l3 = (a1 - a4) ** 2
    l4 = (a2 - a3) ** 2
    l5 = (a2 - a4) ** 2
    l6 = (a3 - a4) ** 2
    v0 = volume_2(l1, l2, l3, l4, l5, l6)
    v1 = volume_2(l1, l2, l4, r1, r2, r3)
    v2 = volume_2(l1, l3, l5, r1, r2, r4)
    v3 = volume_2(l2, l3, l6, r1, r3, r4)
    v4 = volume_2(l4, l5, l6, r2, r3, r4)
    if v0 <= 0: return False
    return v0 - v1 - v2 - v3 - v4 <= 0.0

def check_triangle(a1, a2, a3, r=1.4):
    r1 = (get_radius(a1.element) + r) ** 2
    r2 = (get_radius(a2.element) + r) ** 2
    r3 = (get_radius(a3.element) + r) ** 2
    l1 = (a1 - a2) ** 2
    l2 = (a2 - a3) ** 2
    l3 = (a1 - a3) ** 2
    s0 = area_2(l1, l2, l3)
    s1 = area_2(l2, r2, r3)
    s2 = area_2(l3, r1, r3)
    s3 = area_2(l1, r1, r2)
    if s0 <= 0: return False
    if s1 <= 0: return True
    if s2 <= 0: return True
    if s3 <= 0: return True
    return s0 - s1 - s2 - s3 >= 0.0

def check_triangle2(a1, a2, a3, r=1.4):
    from chempy import cpv
    r1 = (a1.vdw + r) ** 2
    r2 = (a2.vdw + r) ** 2
    r3 = (a3.vdw + r) ** 2
    l1 = cpv.distance(a1.coord, a2.coord) ** 2
    l2 = cpv.distance(a2.coord, a3.coord) ** 2
    l3 = cpv.distance(a1.coord, a3.coord) ** 2
    s0 = area_2(l1, l2, l3)
    s1 = area_2(l2, r2, r3)
    s2 = area_2(l3, r1, r3)
    s3 = area_2(l1, r1, r2)
    if s0 <= 0: return False
    if s1 <= 0: return True
    if s2 <= 0: return True
    if s3 <= 0: return True
    return s0 - s1 - s2 - s3 >= 0.0

if __name__ == "__main__":
    print(area_2(16, 4, 4))
