#!/usr/bin/env python
import numpy as np
# utils for use with numpy and multidimensional arrays

def foo(arr1, arr2, func):
    a = np.vstack((arr1, arr2))
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    return func(b).view(a.dtype).reshape(-1, a.shape[1])

def foo2(arr1, arr2, func):
    b1 = np.ascontiguousarray(arr1).view(np.dtype((np.void, arr1.dtype.itemsize * arr1.shape[1])))
    b2 = np.ascontiguousarray(arr2).view(np.dtype((np.void, arr2.dtype.itemsize * arr2.shape[1])))
    return func(b1, b2).view(arr1.dtype).reshape(-1, arr1.shape[1])

def union(a, b):
    if a == [] or len(a.shape)<2: return b
    if b == [] or len(b.shape)<2: return a
    return foo2(a, b, np.union1d)
def intersect(a, b):
    if a == [] or len(a.shape) < 2: return []
    if b == [] or len(b.shape) < 2: return []
    return foo2(a, b, np.intersect1d)
def diff(a, b):
    if a == [] or len(a.shape) < 2: return []
    if b == [] or len(b.shape) < 2: return a
    return foo2(a, b, np.setdiff1d)


#other
def get_radius(atom_name):
    atom_radii = {"C": 1.7, "N": 1.55, "H": 1.2, "O": 1.52, "S": 1.8}#without hydrogen, simple van der Waals radii
    if (atom_name in atom_radii):
        return atom_radii[atom_name]
    #return 1.5
    #TODO: raise exception if atom is unknown
    #TODO: process situations where hydrogen doesn't appear in pdb file
