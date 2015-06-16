#!/usr/bin/env python

from benchmarkers import *

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
