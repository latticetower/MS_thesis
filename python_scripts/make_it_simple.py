#!/usr/bin/env python

from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

######################################
# decorators (for benchmarking)
######################################
function_timings = {}

def timed(fn):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        result = fn(*args, **kwargs)
        elapsed = timeit.default_timer() - start_time
        logging.debug("Function \"" + fn.__name__ + "\" timed: " + str(elapsed))
        return result
    return wrapper

def accumulated(fn):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        result = fn(*args, **kwargs)
        elapsed = timeit.default_timer() - start_time
        if (function_timings.has_key(fn.__name__)):
            function_timings[fn.__name__][0] += 1
            function_timings[fn.__name__][1] += elapsed
        else:
            function_timings[fn.__name__] = [1, elapsed]
        #logging.debug("Function \"" + fn.__name__ + "\" timed: " + str(elapsed))
        return result
    return wrapper

import atexit


@atexit.register
def print_accumulated():
    logging.debug("Start of accumulated timings")
    for func_name in function_timings.keys():
        [total_calls, total_time] = function_timings[func_name]
        logging.debug(func_name + " : " + str(total_calls) + " calls, elapsed: " + str(total_time) +
            " average: " + str(total_time/total_calls))
    logging.debug("End of accumulated timings")
######################################
# actual functions
######################################

@timed
def read_pdb_info(filename, chain1= 'L', chain2 = 'H'):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    return (
        list(structure[0][chain1].get_atoms()),
        list(structure[0][chain2].get_atoms()) # ,
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['L'].get_atoms())),
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['H'].get_atoms()))
    )

#for given point, construct list of tetrahedra ids
def construct_vertex_dict(triangulation):
    point_no = len(triangulation.points)
    result = [ set() for i in range(point_no)]
    for i in xrange(len(triangulation.simplices)):
        simplex = triangulation.simplices[i]
        for a in simplex:
            result[a].add(i)
    return result
#TODO: in progress
#returns hash of (set -> set)
def construct_edge_dict(triangulation):
    result = {}
    for i in xrange(len(triangulation.simplices)):
        simplex = triangulation.simplices[i]
        for a in get_simplices_2(simplex):
            if not result.has_key(a):
                result[a] = set()
            result[a].add(i)
    return result
def construct_triangle_dict(triangulation):
    result = {}
    for i in xrange(len(triangulation.simplices)):
        simplex = triangulation.simplices[i]
        for a in get_triangles(simplex):
            if not result.has_key(a):
                result[a] = set()
            result[a].add(i)
    return result

# second arg - to test on very-very small subset
# should return triangulation, with some additional objects - neighbours of edge and neighbours of vertex
@timed
def process_chain(points, points_no = 32):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay, KDTree
    p0 = points[: points_no]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    #p = np.array(points[:points_no])
    #print(p[0].get_coord()[1])
    #from scipy.spatial import ConvexHull
    #hull = ConvexHull(p)
    #print(hull.simplices)
    tri = Delaunay(p)
    tree = KDTree(p)
    #print (tri.vertex_neighbor_vertices)
    #print(tri.vertex_to_simplex)
    by_vertex = construct_vertex_dict(tri)
    by_edge = construct_edge_dict(tri)
    by_triangle = construct_triangle_dict(tri)
    #print(tri.vertices)
    return (p0, tri, tree, p, by_vertex, by_edge, by_triangle)

def get_simplices_1(l) : return [np.hstack((l[: i], l[i + 1:])) for i in range(0, len(l))]

def get_simplices_2(l) : return [tuple({l[i], l[j]}) for i in range(0, len(l)) for j in range(len(l)) if i!=j]
def get_triangles(l) : return [tuple({l[i], l[j], l[k]}) for i in range(0, len(l)) for j in range(len(l)) for k in range(0, len(l)) if i!=j and i!=k and j!=k]

#method returns set of convex hull triangles from surface1, whose points lay within given cutoff near surface2
def find_interface_triangles(surface1, surface2, cutoff):
    #print(surface1[1].convex_hull)
    #print(surface1[1].points[surface1[1].convex_hull])
    #print(surface2[2].query(surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff))
    triangles = surface1[1].convex_hull[np.any(surface2[2].query(surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff)[0] < cutoff, axis=1)]
    return triangles
    #print(surface1[1].convex_hull[.count_neighbors(surface1) > 0])
    #surface2[2].count_neighbors(surface1)

#input: array of triangles
#output: aminoacids
def to_aa(triangles, surface):
    if len(triangles) == 0:
        return triangles
    #print surface[0]
    g = np.vectorize(lambda x: surface[0][x].get_parent())
    return np.unique(g(triangles))

import numpy as np
from itertools import combinations, chain
from scipy.misc import comb

def comb_index(n, k):
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)),
                        int, count=count*k)
    return index.reshape(-1, k)

class DTNode:
    def __init__(self, points):
        self.points = points
        self.points_set = set()
        for p in points:
            self.points_set.add(p)
        #print(self.points)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        return isinstance(other, DTNode) and len(np.setdiff1d(self.points, other.points)) == 0
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))
def get_raduis(atom_name):
    atom_radii = {"C": 1.7, "N": 1.55, "H": 1.2, "O": 1.52}#without hydrogen, simple van der Waals radii
    if (atom_name in atom_radii):
        return atom_radii[atom_name]
    #TODO: raise exception if atom is unknown
    #TODO: process situations where hydrogen doesn't appear in pdb file
class DTGraph:
    def __init__(self, surface):
        self.nodes_map = {}
        self.visited = {}
        #print(surface)
        idx = comb_index(4, 3)
        idx2 = comb_index(4, 2)
        check_edge = lambda x1, x2 : get_raduis(x1.element) + get_raduis(x2.element) < x1 - x2
        #print(get_raduis(surface[0][0].element))
        #print(surface[0])
        #print(idx)
        for k in surface[1].simplices:
            tetrahedra_nodes = np.asarray([DTNode(m) for m in k[idx]])
            for [p0, p1] in tetrahedra_nodes[idx2]:
                t_edge = np.vectorize(lambda x: surface[0][x])(p0.intersect(p1))
                if not p0 in self.nodes_map:
                    self.nodes_map[p0] = set()
                if not p0 in self.visited:
                    self.visited[p0] = False
                if not p1 in self.nodes_map:
                    self.nodes_map[p1] = set()
                if not p1 in self.visited:
                    self.visited[p1] = False
                if check_edge(*t_edge):
                    self.nodes_map[p0].add(p1)
                    self.nodes_map[p1].add(p0)

                #print t_edge
        #print(self.nodes_map)
        #print(surface[1].simplices[:, idx])
        #f = np.vectorize(lambda x: set(x.flat))
        #g = np.vectorize(lambda x: x)
        #for k in surface[1].simplices[:, idx]:
        #    print k
    #recursive method for finding all paths in graph
    def get_path(self, start_node):
        self.visited[start_node] = True
        result = []
        for next_node in self.nodes_map[start_node]:
            if not self.visited[next_node]:
                result.append(next_node)
                result.extend(self.get_path(next_node))
        return result
    def find_pockets(self, triangles):
        nodes = set()
        # print(triangles)
        return (np.unique(np.asarray([ x.points for triangle in triangles for x in self.get_path(DTNode(triangle)) ]))[0])


if __name__ == "__main__":
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled simple script")
    pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    cutoff = 32
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    G = DTGraph(surface1)
    nodes = np.setdiff1d(G.find_pockets(triangles), triangles)
    print(to_aa(triangles, surface1))#this returns interface aminoacids with atoms located near surface2
    print(to_aa(nodes, surface1))#this retuns all aminoacids forming pockets except ones shown previously
    print("got surfaces, next should find dist between them")
