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
def process_chain(points, points_no = 9):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay
    p0 = points[: points_no]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    #p = np.array(points[:points_no])
    #print(p[0].get_coord()[1])
    tri = Delaunay(p)
    by_vertex = construct_vertex_dict(tri)
    by_edge = construct_edge_dict(tri)
    by_triangle = construct_triangle_dict(tri)
    #print(tri.vertices)
    return (p0, tri, by_vertex, by_edge, by_triangle)

def get_simplices_1(l) : return [np.hstack((l[: i], l[i + 1:])) for i in range(0, len(l))]

def get_simplices_2(l) : return [tuple({l[i], l[j]}) for i in range(0, len(l)) for j in range(len(l)) if i!=j]
def get_triangles(l) : return [tuple({l[i], l[j], l[k]}) for i in range(0, len(l)) for j in range(len(l)) for k in range(0, len(l)) if i!=j and i!=k and j!=k]

def get_point_projection_to_plane(tri, p):
    v0 = p - tri[0]
    v1 = tri[1] - tri[0]
    v2 = tri[2] - tri[0]
    normal_v1 = vector_to_axis(v1, v0)
    normal_v2 = vector_to_axis(v2, v0)
    projection_v = normal_v1 + normal_v2 - vector_to_axis(normal_v1, normal_v2)
    return projection_v + tri[0]

def select_by_some_value(values,
            select_func = lambda x: min(x, key = lambda y: y[0]),
            compare_func = lambda x, y: x[0] == y[0]
            ):
    some_value = select_func(values)
    some_values = filter(lambda x: compare_func(some_value, x), values)
    return some_values

# returns [(<distance>, tuple(<points from triangle>))]
@accumulated
def find_distance_to_triangle(triangle, surface1, point):
    tri = map(lambda x: surface1[x].get_vector(), triangle)
    p = point.get_vector()
    p0 = get_point_projection_to_plane(tri, p)
    from math import copysign
    line_segments = get_simplices_1(triangle)
    cmp_result = 0
    for line_segment in line_segments:
        v1 = p0 - surface1[line_segment[0]].get_vector()
        v2 = surface1[line_segment[1]].get_vector() - surface1[line_segment[0]].get_vector()
        res = cmp(v1 ** v2, 0)
        if (res == 0):
            norm_v = p - surface1[line_segment[0]] - v1
            return [(norm_v.norm(), tuple(set([line_segment[0], line_segment[1]])))]
        if (cmp_result == 0):
            cmp_result = res
            continue
        if (cmp_result != res):
            values = list(set([
                (   vector_to_axis(
                        surface1[line[1]].get_vector() - surface1[line[0]].get_vector(),
                        p - surface1[line[0]].get_vector()
                    ).norm(),
                    tuple(set([
                        line[0],
                        line[1]
                    ]))
                )
                for line in line_segments
            ]))
            return select_by_some_value(values)
            #return min_value #todo: check if possible more than 2 values with the same dist
        #for all segments cross product has 1 sign - it means that point projection lays inside of triangle
    return [((p0 - p).norm(), tuple(set(triangle)))]

# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
@accumulated
def find_distance_between_line_segments(line_segment1, surface1, line_segment2, surface2):
    p1 = surface1[line_segment1[0]].get_vector()
    v1 = surface1[line_segment1[1]].get_vector() - surface1[line_segment1[0]].get_vector()
    p2 = surface2[line_segment2[0]].get_vector()
    v2 = surface2[line_segment2[1]].get_vector() - surface2[line_segment2[0]].get_vector()
    dist = (v2 * v2) * (v1 * v1) - (v1 * v2) * (v1 * v2)
    if (dist == 0):
        return [(vector_to_axis(v1, p2 - p1).norm(), (tuple(set(line_segment1)), tuple(set(line_segment2))))]
    k2 = ((p2 - p1)*(v1**(v1*v2)-v2**(v1*v1)))/dist
    k1 = ((p2-p1)*v1+(v1*v2)*k2)/(v1*v1)
    vect1 = p1 + v1 ** k1
    vect2 = p2 + v2 ** k2
    return [((vect2 - vect1).norm(), (tuple(set(line_segment1)), tuple(set(line_segment2))))]

# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
@accumulated
def find_distance_between_triangles(triangle1, surface1, triangle2, surface2):
    distances1 = list(set([
            (x[0], ((point), x[1]))
            for point in triangle1
            for x in find_distance_to_triangle(triangle2, surface2, surface1[point])
        ]))
    distances2 = list(set([
            (x[0], (x[1], (point)))
            for point in triangle2
            for x in find_distance_to_triangle(triangle1, surface1, surface2[point])
        ]))
    distances3 = list(set([
                x
                for line_segment1 in get_simplices_1(triangle1)
                for line_segment2 in get_simplices_1(triangle2)
                for x in find_distance_between_line_segments(line_segment1, surface1, line_segment2, surface2)
            ]))
    return list(set(select_by_some_value(distances1 + distances2 + distances3)))

#helper method - returns set of all unvisited neighbours for given tetrahedra
#construct data from tuple returned by process_chain
def get_all_neighbours(tetrahedras, surface, visited, nearest_vertices):
    result = set()
    for nv in nearest_vertices:
        if len(nv) == 3:
            #get all tetrahedras by triangle
            for x in surface[4][nv]: result.add(x)
        if len(nv) == 2:
            #get all tetrahedras by line segment
            for x in surface[3][nv]: result.add(x)
        if len(nv) == 1:
            #get all tetrahedras by point
            for x in surface[2][nv]: result.add(x)
    return result
    #result = set()
    #for t in tetrahedras:
    #    tetrahedra = surface[1].simplices[t]
    #    for v in tetrahedra:
    #        for neighbour in surface[2][v]:
    #            if not visited[neighbour]:
    #                result.add(neighbour)
    #return result

#TODO: in progress
# returns tetrahedras list with minimal length
# returns (distance, tetrahedras)
#tetrahedras is not actual tetrahedras, it's just index in simplices obj
def check_neighbour_ridges(tetrahedras, surface1, given_tetrahedra, surface2, visited_tetrahedras, start_distance_info):
    start_distance = start_distance_info[0][0]
    nearest_vertices = map(lambda x: x[1][0], start_distance_info)
    neighbours = get_all_neighbours(tetrahedras, surface1, visited_tetrahedras, nearest_vertices)
    nearest_tetrahedras = set()
    min_dist = start_distance_info
    for t in neighbours:
        if visited_tetrahedras[t]:
            continue
        neighbour = surface1[1].simplices[t]
        visited_tetrahedras[t] = True
        data = find_distance(neighbour, surface1[0], given_tetrahedra, surface2[0])
        dist = data
        if dist[0][0] < min_dist[0][0]:
            min_dist = dist
            nearest_tetrahedras = set()
        if dist[0][0] == min_dist[0][0]:
            for x in list(dist):
                min_dist.append(x)
            min_dist = list(set(min_dist))
        nearest_tetrahedras.add(t)
    return (min_dist, nearest_tetrahedras)
    #TODO: add list or some object with list of all visited vertices
    #TODO: check if neighbours of all edges are

#TODO: in progress
# method finds tetrahedra obj visible for given point (or tetrahedra?!...)
@accumulated
def find_nearest_ridges(tetrahedras, surface1, given_tetrahedra, surface2):
    start_tetrahedra = tetrahedras[0]
    visited_tetrahedras = [False for x in xrange(len(tetrahedras))]
    visited_tetrahedras[0] = True
    start_distance_info = find_distance(start_tetrahedra, surface1[0], given_tetrahedra, surface2[0])
    start_neighbours = {0}
    while(True):
        (neighbour_distance_info, candidates) = check_neighbour_ridges(start_neighbours, surface1, given_tetrahedra, surface2, visited_tetrahedras, start_distance_info)
        if neighbour_distance_info[0][0] > start_distance_info[0][0]:
            break
        if neighbour_distance_info[0][0] < start_distance_info[0][0]:
            start_distance_info = neighbour_distance_info
            start_neighbours = candidates
        else:
            if (len(candidates) == 0):
                break
            for i in candidates:
                start_neighbours.add(i)
    return (start_distance_info[0][0], start_neighbours)


#TODO: in progress
# method iterates over tetrahedras in triangulation and finds nearest one
# returns (distance, tetrahedra)
def find_some_ridge(tetrahedras, surface1, given_tetrahedra, surface2):
    start_tetrahedra = tetrahedras[0]
    # visited_tetrahedras
    # i'm not sure is it useful or not
    print "ok"

#method returns tuple - nearest distance and corresponding tuple of arrays of atoms from 1st and second tetrahedra
# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
@accumulated
def find_distance(tetrahedra1, surface1, tetrahedra2, surface2):
    dist =  list(set([
            x
            for triangle1 in get_simplices_1(tetrahedra1)
            for triangle2 in get_simplices_1(tetrahedra2)
            for x in find_distance_between_triangles(triangle1, surface1, triangle2, surface2)
        ]))
    ret = select_by_some_value(dist)
    return ret

def get_tetrahedra(simplex, points) : return (map(lambda x: points[x], simplex))

#simple and super-slow: not actual hausdorrf, returns just pairs of nearest atoms now. and converts them to corresponding aa pairs
@timed
def hausdorff_distance(surface1, surface2):
    #distances = [
    #    x
    #    for simplex1 in surface1[1].simplices
    #    for simplex2 in surface2[1].simplices
    #    for x in find_distance(
    #        simplex1, surface1[0],
    #        simplex2, surface2[0]
    #        )
    #]
    distances = [
        x
        for simplex2 in surface2[1].simplices
        for simplex1 in find_nearest_ridges(surface1[1].simplices, surface1, simplex2, surface2)[1]
        for x in find_distance(surface1[1].simplices[simplex1], surface1[0], simplex2, surface2[0])
    ]
    nearest_atoms = select_by_some_value(distances)
    def atoms_to_aa(surface, atoms) :
        return list(set([surface[a].get_parent().get_id()[1] for a in atoms]))
    res = list(set(
        [
            (x[0],(aa1, aa2))
            for x in nearest_atoms
            for aa1 in atoms_to_aa(surface1[0], x[1][0])
            for aa2 in atoms_to_aa(surface2[0], x[1][1])
        ]))
    return res


if __name__ == "__main__":
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled triangulator")
    pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    print("got surfaces, next should find dist between them")
    dist = hausdorff_distance(surface1, surface2)
    print(dist)
