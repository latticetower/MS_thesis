from Bio.PDB import *
import numpy as np
import timeit
import logging

logging.basicConfig(filename='logs/triangulator.log',level=logging.DEBUG)



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

# second arg - to test on very-very small subset
@timed
def process_chain(points, points_no = 9):
    logging.debug("process_chain called for " + str(points_no) + " atoms")
    from scipy.spatial import Delaunay
    p0 = points[: points_no]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    #p = np.array(points[:points_no])
    #print(p[0].get_coord()[1])
    tri = Delaunay(p)
    #print(tri.simplices)
    return (p0, tri.simplices)

def get_simplices_1(l) : return [l[: i] + l[i + 1:] for i in range(0, len(l))]

def get_point_projection_to_plane(tri, p):
    v0 = p - tri[0]
    v1 = tri[1] - tri[0]
    v2 = tri[2] - tri[0]
    normal_v1 = vector_to_axis(v1, v0)
    normal_v2 = vector_to_axis(v2, v0)
    projection_v = normal_v1 + normal_v2 - vector_to_axis(normal_v1, normal_v2)
    return projection_v + tri[0]

def select_by_some_value(l, select_func, reduce_func = lambda x, y: (x[0], (x[1] + y[1]))):
    some_value = select_func(l, key = lambda x: x[0])[0]
    some_values = filter(lambda x: x[0] == some_value, l)
    #return reduce(reduce_func, some_values)
    return some_values

# returns [(<distance>, [<points from triangle>])]
@accumulated
def find_distance_to_triangle(triangle, point):
    tri = map(lambda x: x.get_vector(), triangle)
    p = point.get_vector()
    p0 = get_point_projection_to_plane(tri, p)
    from math import copysign
    line_segments = get_simplices_1(triangle)
    cmp_result = 0
    for line_segment in line_segments:
        v1 = p0 - line_segment[0].get_vector()
        v2 = line_segment[1].get_vector() - line_segment[0].get_vector()
        res = cmp(v1 ** v2, 0)
        if (res == 0):
            norm_v = p - line_segment[0] - v1
            #print((norm_v.norm(), [line_segment[0], line_segment[1]]))
            return [(norm_v.norm(), [line_segment[0], line_segment[1]])]
        if (cmp_result == 0):
            cmp_result = res
            continue
        if (cmp_result != res):
            values = [
                (   vector_to_axis(
                        line[1].get_vector() - line[0].get_vector(),
                        p - line[0].get_vector()
                    ).norm(),
                    [
                        line[0],
                        line[1]
                    ]
                )
                for line in line_segments
            ]
            return select_by_some_value(values, min)
            #return min_value #todo: check if possible more than 2 values with the same dist
        #for all segments cross product has 1 sign - it means that point projection lays inside of triangle
    return [((p0 - p).norm(), triangle)]

# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
@accumulated
def find_distance_between_line_segments(line_segment1, line_segment2):
    p1 = line_segment1[0].get_vector()
    v1 = line_segment1[1].get_vector() - line_segment1[0].get_vector()
    p2 = line_segment2[0].get_vector()
    v2 = line_segment2[1].get_vector() - line_segment2[0].get_vector()
    dist = (v2 * v2) * (v1 * v1) - (v1 * v2) * (v1 * v2)
    if (dist == 0):
        return (vector_to_axis(v1, p2 - p1).norm(), (line_segment1, line_segment2))
    k2 = ((p2 - p1)*(v1**(v1*v2)-v2**(v1*v1)))/dist
    k1 = ((p2-p1)*v1+(v1*v2)*k2)/(v1*v1)
    vect1 = p1 + v1 ** k1
    vect2 = p2 + v2 ** k2
    return [((vect2 - vect1).norm(), (line_segment1, line_segment2))]

# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
@accumulated
def find_distance_between_triangles(triangle1, triangle2):
    distances1 = [
            (x[0], ([point], x[1]))
            for point in triangle1
            for x in find_distance_to_triangle(triangle2, point)
        ]
    distances2 = [
            (x[0], (x[1], [point]))
            for point in triangle1
            for x in find_distance_to_triangle(triangle2, point)
        ]
    distances3 = [
                x
                for line_segment1 in get_simplices_1(triangle1)
                for line_segment2 in get_simplices_1(triangle2)
                for x in find_distance_between_line_segments(line_segment1, line_segment2)
            ]
    return select_by_some_value(distances1 + distances2 + distances3, min)

#method returns tuple - nearest distance and corresponding tuple of arrays of atoms from 1st and second tetrahedra
# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
def find_distance(tetrahedra1, tetrahedra2):
    dist =  [
            x
            for triangle1 in get_simplices_1(tetrahedra1)
            for triangle2 in get_simplices_1(tetrahedra2)
            for x in find_distance_between_triangles(triangle1, triangle2)
        ]
    ret = select_by_some_value(dist, min)
    return ret

def get_tetrahedra(simplex, points) : return (map(lambda x: points[x], simplex))

#simple and super-slow: not actual hausdorrf, returns just pairs of nearest atoms now. and converts them to corresponding aa pairs
@timed
def hausdorff_distance(surface1, surface2):
    distances = [
        x
        for simplex1 in surface1[1]
        for simplex2 in surface2[1]
        for x in find_distance(
            get_tetrahedra(simplex1, surface1[0]),
            get_tetrahedra(simplex2, surface2[0])
            )
    ]
    nearest_atoms = select_by_some_value(distances, min)
    def atoms_to_aa(atoms) :
        return list(set([a.get_parent().get_id()[1] for a in atoms]))
    res = list(set(
        [
            (x[0],(aa1, aa2))
            for x in nearest_atoms
            for aa1 in atoms_to_aa(x[1][0])
            for aa2 in atoms_to_aa(x[1][1])
        ]))
    #print res
    return res


if __name__ == "__main__":
    logging.debug("============\nCalled triangulator")
    pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    print("got surfaces, next should find dist between them")
    dist = hausdorff_distance(surface1, surface2)
    print(dist)
