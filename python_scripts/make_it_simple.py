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

def group_elements(data):
    coils = ['-', 'B', 'T', 'E']
    last_coil = None
    start_coil_pos = -1
    protein_coils = [] #contains tuples - (start position, end upper bound position, coil type)
    last_key = -1
    for (key, value) in data:
        if value in coils:
            if last_coil != value:
                if last_coil != None:
                    protein_coils.append((start_coil_pos, last_key, value))
                last_coil = value
                start_coil_pos = key
        else:
            if last_coil != None:
                protein_coils.append((start_coil_pos, last_key, last_coil))
            last_coil = None
        last_key = key
    if last_coil != None:
        protein_coils.append((start_coil_pos, last_key, value))
    return protein_coils

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

@timed
def read_dssp_info(filename,
                    chain1= 'L', chain2 = 'H',
                    key_transf1 = lambda x: long(x.get_id()[1]),
                    key_transf2 = lambda x: x.get_parent().get_id()):
    chains_data = {}
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    model = structure[0]
    dssp = DSSP(model, filename, dssp='mkdssp')
    # DSSP data is accessed by a tuple (chain_id, res_id)
    #a_key = list(dssp)[2]
    #last_key = -1
    #for (key, value, v, e, r,u) in dssp:
    #    print key.get_parent().get_id()
    #return
    for chain in (chain1, chain2):
        data = [
                (key_transf1(key), value)
                for (key, value, v, e, r,u) in dssp
                if key_transf2(key) == chain
            ]
        chains_data[chain] = group_elements(data)
    return chains_data


# second arg - to test on very-very small subset
# should return triangulation, with some additional objects - neighbours of edge and neighbours of vertex
@timed
def process_chain(points, points_no = 0):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    print("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay, KDTree
    p0 = points[: (points_no if points_no > 0 else len(points))]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    tri = Delaunay(p)
    tree = KDTree(p)
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
    triangles = surface1[1].convex_hull[np.any(surface2[2].query(
        surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff
        )[0] < cutoff, axis = 1)]
    return triangles

#input: array of triangles
#output: aminoacids
def to_aa(triangles, surface):
    if len(triangles) == 0:
        return triangles
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
                        int, count = count * k)
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
        assert(len(self.points) == len(other.points) == 3)
        cond3 = len(np.setdiff1d(self.points, other.points)) == 0
        return isinstance(other, DTNode) and isinstance(self, DTNode) and cond3
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))

def get_raduis(atom_name):
    atom_radii = {"C": 1.7, "N": 1.55, "H": 1.2, "O": 1.52, "S": 1.8}#without hydrogen, simple van der Waals radii
    if (atom_name in atom_radii):
        return atom_radii[atom_name]
    #return 1.5
    #TODO: raise exception if atom is unknown
    #TODO: process situations where hydrogen doesn't appear in pdb file
class DTGraph:
    def __init__(self, surface):
        self.surface = surface
        self.path_result = False
        self.nodes_map = {}
        self.visited = {}
        self.nearest_triangles = {}
        idx = comb_index(4, 3)
        idx2 = comb_index(4, 2)
        self.ch_nodes = set([DTNode(t) for t in surface[1].convex_hull])
        for node in self.ch_nodes:
            self.nearest_triangles[node] = node
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
        for s in self.ch_nodes:
            self.visited[s] = True #to avoid selection of all outer non-convex area #won't help
    def is_ch_simplex(self, node):
        assert(isinstance(node, DTNode))
        return node in self.ch_nodes
    def dist(self, node1, node2, distance_func):
        #print(node1)
        l1 = np.vectorize(lambda x: self.surface[0][x])(node1.points)
        l2 = np.vectorize(lambda x: self.surface[0][x])(node2.points)
        #d0 = ([self.surface[0][p1] for p1 in node1.points for p2 in node2.points]).pop()
        #sprint d0
        #print d0
        d1 = min([distance_func(p1, p2) for p1 in l1 for p2 in l2])
        #check_edge = lambda x1, x2 : x1.vdw + x2.vdw < cpv.distance(x1.coord, x2.coord)
        #print d1
        return d1
    def find_nearest_node(self, next_node, distance_func):
        if next_node in self.nearest_triangles:
            return self.nearest_triangles[next_node]
        self.nearest_triangles[next_node] = next(iter(self.ch_nodes))
        for node in self.ch_nodes:
            #print node
            if self.dist(next_node, node, distance_func) < self.dist(next_node, self.nearest_triangles[next_node], distance_func):
                self.nearest_triangles[next_node] = node
    def check_dist(self, triangle, start_nodes):
        #print(self.nearest_triangles)
        #print(start_nodes)
        return self.nearest_triangles[triangle] in start_nodes
    # method adds all child nodes and returns current node (from queue)
    def get_path_fragment(self, check_edge, start_nodes, distance_func):
        if len(self.queue) == 0:
            return []
        while True:
            start_node = self.queue.pop(0)
            if not self.visited[start_node] or self.is_ch_simplex(start_node):
                break
            if (len(self.queue)==0):
                return []
        self.visited[start_node] = True
        if not start_node in self.nodes_map:
            return []
        result = set() #start_node
        for next_node in self.nodes_map[start_node]:
            edge = self.nodes_map[start_node][next_node]
            if next_node in self.visited:
                self.find_nearest_node(next_node, distance_func)
                if check_edge(*edge) and not self.is_ch_simplex(next_node) and self.check_dist(next_node, start_nodes):
                    if not self.visited[next_node]:
                        if not next_node in self.queue:
                            self.queue.append(next_node)
        result.add(start_node)
        return result
    def find_pockets(self, triangles, check_edge, distance_func):
        start_nodes = set([DTNode(triangle) for triangle in triangles])
        self.queue = [DTNode(triangle) for triangle in triangles]
        data = set()
        while (len(self.queue) > 0):
            #print("".join(["T" if self.visited[DTNode(t)] else "F" for t in triangles]))
            #print("".join(["T" if self.visited[t] else "F" for t in self.queue]))
            for x in self.get_path_fragment(check_edge, start_nodes, distance_func):
                data.add(x)
        result = np.asarray([x.points for x in data])
        if(len(result) == 0):
            return result
        return (result)
    # this version of pocket finder should use limited Dijkstra's search:
    #
    #def find_pockets2(self, triangles, check_edge):

class CHNode:
    def __init__(self, points):
        self.points = points
        self.points_set = set()
        for p in points:
            self.points_set.add(p)
        assert(len(points) == 2)
        assert(len(self.points_set) == 2)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        assert(len(self.points) == len(other.points) == 2)
        #try:
        cond2 = len(np.setdiff1d(self.points, other.points)) == 0
        return isinstance(other, CHNode) and isinstance(self, CHNode) and cond2
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))

class CHGraph:
    def __init__(self, surface):
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

#following method adds coil aminoacids to aminoacids touched by interface_triangles atoms
def extend_to_coils(interface_triangles, chain_info, surface, get_res_seq = lambda x: x.get_parent().get_id()[1]):
    if len(interface_triangles) == 0:
        return interface_triangles
    g = np.vectorize(lambda x: long(get_res_seq(surface[0][x])))
    result =  np.unique(g(interface_triangles))
    coil_fragments = set()
    distinct_aa = set()
    for aa_id in result:
        coil = [c for c in chain_info if aa_id >= c[0] and aa_id <= c[1] ]
        if len(coil) > 0:
            for c in coil:
                coil_fragments.add(c)
        else:
            distinct_aa.add(aa_id)
    for fragment in coil_fragments:
        for k in xrange(fragment[0], fragment[1]):
            distinct_aa.add(k)
    return np.unique(np.asarray(list(distinct_aa)))

if __name__ == "__main__":
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled simple script")
    pair_of_chains = read_pdb_info(os.path.abspath('../test_data/2OSL.pdb'))
    chains_ss_info = read_dssp_info(os.path.abspath('../test_data/2OSL.pdb'))
    #print chains_ss_info
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    cutoff = 5.0
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
    def distance_func(x1, x2):
        return x1 - x2
    #check_edge = lambda x1, x2 : (get_raduis(x1.element) + get_raduis(x2.element) < x1 - x2)
    #check_edge2 = lambda x1, x2 : True
    G = DTGraph(surface1)
    nodes = np.setdiff1d(G.find_pockets(triangles, check_edge, distance_func), triangles)
    #nodes = G.find_pockets(triangles, check_edge)
    print(nodes)
    aa_with_coils = extend_to_coils(triangles, chains_ss_info['L'], surface1)
    print aa_with_coils
    #    if (len(nodes) > 0):
    #print(to_aa(nodes, surface1))#this retuns all aminoacids forming pockets except ones shown previously
    print("end of processing")
    #C = CHGraph(surface1)
