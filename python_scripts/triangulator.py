from Bio.PDB import *
import numpy as np

def read_pdb_info(filename):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    return (
        list(map(lambda x : (x.get_vector()._ar), structure[0]['L'].get_atoms())),
        list(map(lambda x : (x.get_vector()._ar), structure[0]['H'].get_atoms()))
    )
def process(points):
    from scipy.spatial import Delaunay
    points_no = 31
    p = np.array(points[:points_no])
    tri = Delaunay(p)
    return (p, tri.simplices)

def find_distace_between_points(p1, p2):
    sqrt(sum(map(lambda x : x[0] * x[1], zip(p1, p2))))
def find_distace(tetrahedra1, tetrahedra2):
    print("ok")

def hausdorff_distance(surface1, surface2):
    for simplex1 in surface1[1]:
        for simplex2 in surface2[1]:
            find_distace(surface1[0][simplex1], surface2[0][simplex2])
    print("ok")

pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')

surface1 = process(pair_of_chains[0])
surface2 = process(pair_of_chains[1])

dist = hausdorff_distance(surface1, surface2)
print(dist)
