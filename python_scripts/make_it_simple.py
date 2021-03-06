#!/usr/bin/env python

from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

from benchmarkers import *
from utils import *
######################################
# actual functions
######################################

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

def group_all_elements(data):
    grouped_data = {}
    for (pos, key) in data:
        if not key in grouped_data:
            grouped_data[key] = set()
        grouped_data[key].add(pos)
    return grouped_data

#todo: add water molecules filtering
@timed
def read_pdb_info(filename, chain1= 'L', chain2 = 'H'):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    if chain2 == None:
        chain_names = map(lambda x: x.id, structure[0])
        if len(chain_names) != 2:
            print "there are more than 2 chains, can't decide"
            return None
        if chain1 != None:
            chain2 = chain_names[1] if chain1 == chain_names[0] else chain_names[0]
    chain1_data = [atom for atom in structure[0][chain1].get_atoms() if atom.get_parent().get_id()[0] != 'W']
    chain2_data = [atom for atom in structure[0][chain2].get_atoms() if atom.get_parent().get_id()[0] != 'W']
    chain1_residues = {atom.get_parent().get_id()[1] for atom in chain1_data}
    #print chain1_residues
    return (
        chain1_data,
        chain2_data,
        len(chain1_residues) # ,
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['L'].get_atoms())),
        #list(map(lambda x : (x.get_vector()._ar), structure[0]['H'].get_atoms()))
    )

@timed
def read_pdb_chains_info(filename):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    return [ chain.get_id() for chain in structure[0] ]

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
    print p
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

#method returns set of all atoms from surface1, whose points lay within given cutoff near surface2
def find_all_by_cutoff(surface1, surface2, cutoff):
    #TODO: check and write it
    triangles = surface1[1].convex_hull[np.any(surface2[2].query(
        surface1[1].points, 1, 0, 2, cutoff
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
    @accumulated
    def __init__(self, points):
        self.points = points
        self.points_set = set()
        for p in points:
            self.points_set.add(p)
        #assert(len(points)==3)
        #assert(len(self.points_set)==3)
        #print(self.points)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        #assert(len(self.points) == len(other.points) == 3)
        cond3 = len(np.setdiff1d(self.points, other.points)) == 0
        return isinstance(other, DTNode) and isinstance(self, DTNode) and cond3
    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        return hash(tuple(self.points_set))

def get_radius(atom_name):
    atom_radii = {"C": 1.7, "N": 1.55, "H": 1.2, "O": 1.52, "S": 1.8}#without hydrogen, simple van der Waals radii
    if (atom_name in atom_radii):
        return atom_radii[atom_name]
    #return 1.5
    #TODO: raise exception if atom is unknown
    #TODO: process situations where hydrogen doesn't appear in pdb file

class DTGraph:
    @timed
    def __init__(self, surface):
        self.surface = surface
        self.path_result = False
        self.nodes_map = {}
        self.visited = set()
        self.nearest_triangles = {}
        self.idx = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
        #comb_index(4, 3)
        #print idx
        #idx2 = comb_index(4, 2)
        self.idx2 = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
        #idx2_3 = comb_index(3, 2)
        self.idx2_3 = np.array([[0, 1], [0, 2], [1, 2]])
        self.ch_neighbours = {}
        #print(surface[1].convex_hull)
        self.ch_nodes = set()
        for t in surface[1].convex_hull:
            node = DTNode(t)
            self.ch_nodes.add(node)
            for x in t[self.idx2_3]:
                if tuple(set(x)) not in self.ch_neighbours:
                    self.ch_neighbours[tuple(set(x))] = []
                self.ch_neighbours[tuple(set(x))].append(node)
        for node in self.ch_nodes:
            self.nearest_triangles[node] = node
        """
        for tetrahedra in surface[1].simplices:
            assert(len(tetrahedra) == 4)
            tetrahedra_nodes = np.asarray([DTNode(m) for m in tetrahedra[self.idx]])
            for [p0, p1] in tetrahedra_nodes[self.idx2]:
                assert(p0 != p1)
                t_edge = np.vectorize(lambda x: surface[0][x])(p0.intersect(p1))
                assert(len(t_edge) == 2)
                #if not p0 in self.visited:
                self.visited[p0] = False
                #if not p1 in self.visited:
                self.visited[p1] = False
                #if check_edge(*t_edge):
                if not p0 in self.nodes_map:
                    self.nodes_map[p0] = {}
                if not p1 in self.nodes_map:
                    self.nodes_map[p1] = {}
                self.nodes_map[p0][p1] = t_edge
                self.nodes_map[p1][p0] = t_edge
        """
        self.visited.update(self.ch_nodes) #to avoid selection of all outer non-convex area #won't help
    #helper method to build graph and init data
    @timed
    def build_graph(self, check_edge, distance_func, check_simplex):
        self.nodes_map = {}
        self.edges_checker = {}
        self.simplex_checker = {}
        for tetrahedra in self.surface[1].simplices:
            #assert(len(k) == 4)
            tetrahedra_nodes = np.asarray([DTNode(m) for m in tetrahedra[self.idx]])
            for [node0, node1] in tetrahedra_nodes[self.idx2]:
                t_edge = tuple(set(node0.intersect(node1)))
                if not node0 in self.nodes_map:
                    self.nodes_map[node0] = {}
                if not node1 in self.nodes_map:
                    self.nodes_map[node1] = {}
                if not t_edge in self.edges_checker:
                    self.edges_checker[t_edge] = check_edge(*t_edge)
                if not tuple(node0.points_set) in self.simplex_checker:
                    self.simplex_checker[tuple(node0.points_set)] = check_simplex(tuple(node0.points_set))
                if not tuple(node1.points_set) in self.simplex_checker:
                    self.simplex_checker[tuple(node1.points_set)] = check_simplex(tuple(node1.points_set))
                #condition = self.simplex_checker[tuple(node0.points_set)] #self.edges_checker[t_edge]
                if self.simplex_checker[tuple(node0.points_set)]:
                    self.nodes_map[node0][node1] = t_edge
                if self.simplex_checker[tuple(node1.points_set)]:
                    self.nodes_map[node1][node0] = t_edge
                    #self.nodes_map[node0][node1] = t_edge
                    #self.nodes_map[node1][node0] = t_edge
        #print 1
    def is_ch_simplex(self, node):
        #assert(isinstance(node, DTNode))
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
        #print (start_nodes)
        #print neighbour_nodes
        if next_node in self.nearest_triangles:
            return self.nearest_triangles[next_node]
        self.nearest_triangles[next_node] = next(iter(self.ch_nodes))
        for node in self.start_nodes:
            if self.dist(next_node, node, distance_func) < self.dist(next_node, self.nearest_triangles[next_node], distance_func):
                self.nearest_triangles[next_node] = node
        for node in self.neighbour_nodes:
            if self.dist(next_node, node, distance_func) < self.dist(next_node, self.nearest_triangles[next_node], distance_func):
                self.nearest_triangles[next_node] = node
    def check_dist(self, triangle):
        #print(self.nearest_triangles)
        #print(start_nodes)
        return self.nearest_triangles[triangle] in self.start_nodes
    # method adds all child nodes and returns current node (from queue)
    def get_path_fragment(self, check_edge, distance_func):
        if len(self.queue) == 0:
            return []
        while True:
            start_node = self.queue.pop()
            if not start_node in self.visited or self.is_ch_simplex(start_node):
                break
            if (len(self.queue)==0):
                return []
        self.visited.add(start_node)
        if not start_node in self.nodes_map:
            return []
        #result = set() #start_node
        nodes = set()
        for next_node in self.nodes_map[start_node]:
            #edge = self.nodes_map[start_node][next_node]
            self.find_nearest_node(next_node, distance_func)
            #check_edge(*edge) and
            if not self.is_ch_simplex(next_node) and self.check_dist(next_node):
                if not next_node in self.visited:
                    #if not next_node in self.queue:
                    nodes.add(next_node)
        #for node in nodes:
        #if not next_node in self.queue:
        self.queue.update(nodes)
        #result.add(start_node)
        return start_node
    def get_neighbours(self, start_nodes):
        result = np.unique(np.asarray([k for x in start_nodes for pair in x.points[self.idx2_3] for k in self.ch_neighbours[tuple(set(pair))] if k != x and not k in start_nodes]))
        #print len(result)
        #print len(start_nodes)
        #result = {x for x in result}
        #print len(result)
        return result
    @timed
    def find_pockets(self, triangles, check_edge, distance_func, check_simplex):
        self.start_nodes = {DTNode(triangle) for triangle in triangles}
        self.neighbour_nodes = self.get_neighbours(self.start_nodes)
        self.build_graph(check_edge, distance_func, check_simplex)
        #print(len(self.start_nodes))
        #print(len(self.neighbour_nodes))
        self.queue = {DTNode(triangle) for triangle in triangles}
        data = set()
        while (len(self.queue) > 0):
            #print("".join(["T" if self.visited[DTNode(t)] else "F" for t in triangles]))
            #print("".join(["T" if self.visited[t] else "F" for t in self.queue]))
            x = self.get_path_fragment(check_edge, distance_func)
            data.add(x)
        result = np.asarray([x.points for x in data])
        if(len(result) == 0):
            return result
        return result
    # this version of pocket finder should use limited Dijkstra's search:
    #
    #def find_pockets2(self, triangles, check_edge):

class CHNode:
    def __init__(self, points):
        self.points = points
        self.points_set = set()
        for p in points:
            self.points_set.add(p)
        #assert(len(points) == 2)
        #assert(len(self.points_set) == 2)
    def intersect(self, other):
        return np.intersect1d(self.points, other.points)
    def __eq__(self, other):
        #assert(len(self.points) == len(other.points) == 2)
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
def extend_to_coils(interface_triangles,
                    chain_info,
                    surface,
                    get_res_seq = lambda x: x.get_parent().get_id()[1]
                    ):
    #logging.debug("coils extension started")
    if len(interface_triangles) == 0:
        return interface_triangles
    g = np.vectorize(lambda x: int(get_res_seq(surface[0][x])))
    result =  np.unique(g(interface_triangles))
    coil_fragments = set()
    distinct_aa = set()
    #logging.debug("print chains information")
    #print chain_info
    for aa_id in result:
        #print "aa_id "
        #print aa_id
        coil = [c for c in chain_info if aa_id >= c[0] and aa_id <= c[1] ]
        #logging.debug("aa {0}: {1}".format(aa_id, coil))
        #print coil
        if len(coil) > 0:
            for c in coil:
                coil_fragments.add(c)
        #else:
        #distinct_aa.add(aa_id)
    for fragment in coil_fragments:
        for k in range(fragment[0], fragment[1]+1):
            distinct_aa.add(k)
    return np.unique(np.asarray(list(distinct_aa)))

@timed
def main_func(pdb_filename, chain1, chain2, cutoff = 5.0, sas_radius=1.4):
    pair_of_chains = read_pdb_info(pdb_filename, chain1, chain2)
    if pair_of_chains == None:
        return None
    chains_ss_info = read_dssp_info(pdb_filename, chain1, chain2)
    if chains_ss_info == None:
        return None
    #print chains_ss_info
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    #print(to_aa(triangles, surface1))
    triangles = extend_interface_1(triangles, surface1)
    #print(to_aa(triangles, surface1))#this returns interface aminoacids with atoms located near surface2
    # custom function: returns True if can travel from one triangle to another
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
    #check_edge = lambda x1, x2 : (get_radius(x1.element) + get_radius(x2.element) < x1 - x2)
    #check_edge2 = lambda x1, x2 : True
    G = DTGraph(surface1)
    pockets = G.find_pockets(triangles, check_edge, distance_func, check_simplex)
    print "pockets"
    print pockets
    nodes = diff(pockets, triangles)
    print "nodes"
    print nodes
    #nodes = G.find_pockets(triangles, check_edge)
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

def main_func1(pdb_filename, chain1, chain2, cutoff = 5.0):
    pair_of_chains = read_pdb_info(pdb_filename, chain1, chain2)
    if pair_of_chains == None:
        return None
    chains_ss_info = read_dssp_info(pdb_filename, chain1, chain2)
    if chains_ss_info == None:
        return None
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    chain_length = 0
    return (to_aa(triangles, surface1), #this retuns all aminoacids forming pockets except ones shown previously
        chain_length
    )


def find_by_cutoff(pdb_filename, chain1, chain2, cutoff=5.0):
    pair_of_chains = read_pdb_info(pdb_filename, chain1, chain2)
    if pair_of_chains == None:
        return None
    #chains_ss_info = read_dssp_info(pdb_filename, chain1, chain2)
    #if chains_ss_info == None:
    #    return None
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    pairs = surface1[2].query_ball_tree(surface2[2], cutoff)
    chain_length = pair_of_chains[2]
    #print pairs
    atoms = [surface1[0][i] for i in range(0,len(surface1[2].data)) if len(pairs[i]) > 0]
    return (np.unique(np.asarray(list({atom.get_parent() for atom in atoms}))),
        chain_length
        )
    #triangles = find_all_by_cutoff(surface1, surface2, cutoff)
    #def check_edge(x1, x2):
    #    dist = (x1 - x2) - (get_radius(x1.element) + get_radius(x2.element))
    #    return dist > 0
    #def distance_func(x1, x2):
    #    return x1 - x2
    #return to_aa(triangles, surface1)
    #this retuns all aminoacids forming pockets except ones shown previously


if __name__ == "__main__":
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled simple script")
    res = main_func(os.path.abspath('../test_data/2OSL.pdb'), 'L', 'H')
    print res[0]
    print("end of processing")
