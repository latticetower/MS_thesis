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
    #model = structure[0]
    #dssp = DSSP(model, filename, dssp='mkdssp')
    # DSSP data is accessed by a tuple (chain_id, res_id)
    #a_key = list(dssp)[2]
    #for (key, value, v, e, r,u) in dssp:
    #    print value
    # residue object, secondary structure, solvent accessibility,
    # relative accessiblity, phi, psi
    return (
        list(structure[0][chain1].get_atoms()),
        list(structure[0][chain2].get_atoms()) # ,
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['L'].get_atoms())),
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['H'].get_atoms()))
    )


# second arg - to test on very-very small subset
# should return triangulation, with some additional objects - neighbours of edge and neighbours of vertex
@timed
def process_chain(points, points_no = 0):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    print("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay, KDTree
    p0 = points[: (points_no if points_no > 0 else len(points))]
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
    #print(tri.vertices)
    return (p0, tri, tree, p)

#method returns chain with atom info
def process_atom_chain(points, points_no = 0):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay, KDTree
    p0 = points[: (points_no if points_no > 0 else len(points))]
    p = np.array(map(lambda x : (x.coord), p0))
    tri = Delaunay(p)
    tree = KDTree(p)
    return (p0, tri, tree, p)

def get_simplices_1(l) : return [np.hstack((l[: i], l[i + 1:])) for i in range(0, len(l))]

def get_simplices_2(l) : return [tuple({l[i], l[j]}) for i in range(0, len(l)) for j in range(len(l)) if i!=j]
def get_triangles(l) : return [tuple({l[i], l[j], l[k]}) for i in range(0, len(l)) for j in range(len(l)) for k in range(0, len(l)) if i!=j and i!=k and j!=k]

#method returns set of convex hull triangles from surface1, whose points lay within given cutoff near surface2
def find_interface_triangles(surface1, surface2, cutoff):
    #print(surface1[1].convex_hull)
    #print(surface1[1].points[surface1[1].convex_hull])
    #print(np.any(surface2[2].query(surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff)[0]<cutoff, axis = 1))
    triangles = surface1[1].convex_hull[np.any(surface2[2].query(
        surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff
        )[0] < cutoff, axis=1)]
    #print(triangles)
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

def to_aa2(triangles, surface):
    #print(triangles)
    if len(triangles) == 0:
        return []
    g = np.vectorize(lambda x: surface[0][x])
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
        assert(len(points)==3)
        assert(len(self.points_set)==3)
        #print(self.points)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        assert(len(self.points) == len(other.points) ==3)
        #try:
        cond3 = len(np.setdiff1d(self.points, other.points)) == 0
        return isinstance(other, DTNode) and isinstance(self, DTNode) and cond3
        #except Exception:
        #    print(self.points)
        #    print(other.points)
        #    return True
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))

def get_raduis(atom_name):
    atom_radii = {"C": 1.7, "N": 1.55, "H": 1.2, "O": 1.52, "S": 1.8}#without hydrogen, simple van der Waals radii
    if (atom_name in atom_radii):
        return atom_radii[atom_name]
    #print atom_name
    #return 1.5
    #TODO: raise exception if atom is unknown
    #TODO: process situations where hydrogen doesn't appear in pdb file
class DTGraph:
    def __init__(self, surface):
        self.surface = surface
        self.path_result = False
        self.nodes_map = {}
        self.visited = {}
        #print(surface)
        idx = comb_index(4, 3)
        idx2 = comb_index(4, 2)
        #print(get_raduis(surface[0][0].element))
        #print(surface[0])
        #print(idx)
        self.ch_nodes = set([DTNode(t) for t in surface[1].convex_hull])
        for k in surface[1].simplices:
            assert(len(k) == 4)
            tetrahedra_nodes = np.asarray([DTNode(m) for m in k[idx]])
            for [p0, p1] in tetrahedra_nodes[idx2]:
                assert(p0 != p1)
                t_edge = np.vectorize(lambda x: surface[0][x])(p0.intersect(p1))
                assert(len(t_edge) == 2)
                if not p0 in self.visited:
                    self.visited[p0] = False
                if not p1 in self.visited:
                    self.visited[p1] = False
                #if check_edge(*t_edge):
                if not p0 in self.nodes_map:
                    self.nodes_map[p0] = {}
                if not p1 in self.nodes_map:
                    self.nodes_map[p1] = {}
                self.nodes_map[p0][p1] = t_edge
                self.nodes_map[p1][p0] = t_edge
                #print t_edge
        #print(self.nodes_map)
        #print([k for k in self.visited if not isinstance(k, DTNode)])
        #print(surface[1].simplices[:, idx])
        #f = np.vectorize(lambda x: set(x.flat))
        #g = np.vectorize(lambda x: x)
        #for k in surface[1].simplices[:, idx]:
        #    print k
    def is_ch_simplex(self, node):
        assert(isinstance(node, DTNode))
        return node in self.ch_nodes
    # method adds all child nodes and returns current node (from queue)
    def get_path_fragment(self, check_edge):
        if len(self.queue) == 0:
            return []
        while True:
            start_node = self.queue.pop()
            if not self.visited[start_node]:
                break
            if (len(self.queue)==0):
                return []
        self.visited[start_node] = True
        if not start_node in self.nodes_map:
            return []
        result = set() #start_node
        for next_node in self.nodes_map[start_node]:
            #print(next_node)
            edge = self.nodes_map[start_node][next_node]
            #print(check_edge(*edge))
            #print(self.visited[next_node])
            #print(next_node in self.visited)
            if next_node in self.visited:
                #print(1)
                if check_edge(*edge) and not self.is_ch_simplex(next_node):
                    if not self.visited[next_node]:
                        #print ("+add new node to queue")
                        #print(next_node)
                        if not next_node in self.queue:
                            self.queue.add(next_node)
                #result.add(next_node)
        result.add(start_node)
        return result
    def find_pockets(self, triangles, check_edge):
        #print(triangles)
        self.queue = set([DTNode(triangle) for triangle in triangles])
        data = set()
        while (len(self.queue) > 0):
            #print("".join(["T" if self.visited[DTNode(t)] else "F" for t in triangles]))
            #print("".join(["T" if self.visited[t] else "F" for t in self.queue]))
            for x in self.get_path_fragment(check_edge):
                data.add(x)
        result = np.asarray([x.points for x in data])
        if(len(result) == 0):
            return result
        return (result)
        #return np.unique(result)
        #return #(np.unique(
        #return np.asarray([ x.points for triangle in triangles for x in self.get_path(DTNode(triangle)) ])
        #)[0])
class CHNode:
    def __init__(self, points):
        self.points = points
        self.points_set = set()
        for p in points:
            self.points_set.add(p)
        assert(len(points)==2)
        assert(len(self.points_set)==2)
        #print(self.points)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        assert(len(self.points) == len(other.points) ==2)
        #try:
        cond2 = len(np.setdiff1d(self.points, other.points)) == 0
        return isinstance(other, CHNode) and isinstance(self, CHNode) and cond2
        #except Exception:
        #    print(self.points)
        #    print(other.points)
        #    return True
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))
class CHGraph:
    def __init__(self, surface):
        #print(surface1[1].convex_hull)
        self.edges_list = {}
        self.triangles_list = {}
        idx = comb_index(3, 2)
        for triangle in surface1[1].convex_hull:
            for [l1, l2, l3] in triangle[idx]:
                line = CHNode(l1)
                if not line in self.edges_list:
                    self.edges_list[line] = set()
                self.edges_list[line].add(l2)
                self.edges_list[line].add(l3)
                line = CHNode(l2)
                if not line in self.edges_list:
                    self.edges_list[line] = set()
                self.edges_list[line].add(l1)
                self.edges_list[line].add(l3)
                line = CHNode(l3)
                if not line in self.edges_list:
                    self.edges_list[line] = set()
                self.edges_list[line].add(l1)
                self.edges_list[line].add(l2)
                #self.edges_list[line].add(DTNode(triangle))
        #for edge in self.edges_list:
        #    for self.edges_list[edge]

def extend_interface_1(triangles, surface):
    #print 1
    #aa = [atom for t in triangles for x in t for atom in surface[0][x].get_parent()]
    #print aa
    #tri = triangles
    #pp = surface[1].vertex_neighbor_vertices
    #ppp = np.unique(np.asarray([
    #    point
    #    for p in np.unique(triangles)
    #    for point in pp[1][pp[0][p] : pp[0][p + 1]]
    #    ]))
    #print np.intersect1d(ppp, triangles)
    return triangles
if __name__ == "__main__":
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled simple script")
    pair_of_chains = read_pdb_info(os.path.abspath('../test_data/2OSL.pdb'))
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    cutoff = 5.0
    cutoff = 33
    #while cutoff < 42.0:
    #    print(cutoff)
    #    cutoff += 1
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    #print(to_aa(triangles, surface1))
    triangles = extend_interface_1(triangles, surface1)
    print(to_aa(triangles, surface1))#this returns interface aminoacids with atoms located near surface2
    def check_edge(x1, x2):
        dist = (x1 - x2) - (get_raduis(x1.element) + get_raduis(x2.element))
        #print(dist)
        return dist > 0
    #check_edge = lambda x1, x2 : (get_raduis(x1.element) + get_raduis(x2.element) < x1 - x2)
    #check_edge2 = lambda x1, x2 : True
    G = DTGraph(surface1)
    nodes = np.setdiff1d(G.find_pockets(triangles, check_edge), triangles)
    #nodes = G.find_pockets(triangles, check_edge)
    print(nodes)
    #    if (len(nodes) > 0):
    #print(to_aa(nodes, surface1))#this retuns all aminoacids forming pockets except ones shown previously
    print("end of processing")
    #C = CHGraph(surface1)
