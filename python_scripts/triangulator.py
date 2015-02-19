from Bio.PDB import *
import numpy as np

def read_pdb_info(filename):
    parser = PDBParser()
    structure = parser.get_structure('', filename)
    #print(list(structure[0]['L'].get_residues())[45].get_id()[1])
    #print(list(structure[0]['L'].get_residues())[45].get_id()[1])
    return (
        list(structure[0]['L'].get_atoms()),
        list(structure[0]['H'].get_atoms()),
        list(map(lambda x : (x.get_vector()._ar), structure[0]['L'].get_atoms())),
        list(map(lambda x : (x.get_vector()._ar), structure[0]['H'].get_atoms()))
    )
# second arg - to test on very-very small subset
def process_chain(points, points_no = 9):
    from scipy.spatial import Delaunay
    p0 = points[: points_no]
    p = np.array(map(lambda x : (x.get_vector()._ar), p0))
    #p = np.array(points[:points_no])
    #print(p[0].get_coord()[1])
    tri = Delaunay(p)
    #print(tri.simplices)
    return (p0, tri.simplices)

def get_simplices_1(l) : return [l[: i] + l[i + 1:] for i in range(0, len(l))]

def find_distance_between_points(p1, p2):
    sqrt(sum(map(lambda x : x[0] * x[1], zip(p1, p2))))

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
            print((norm_v.norm(), [line_segment[0], line_segment[1]]))
            return (norm_v.norm(), [line_segment[0], line_segment[1]])
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
        #for all segments cross product has 1 sign - it means that point projection lies inside of triangle
    return [(p0 - p).norm(), triangle]
    #from math import copysign
    #sign1 = (copysign(1, (v1 ** projection).norm()))
    #sign2 = (copysign(1, (projection**v2).norm()))
    #print(sign1 == sign2)
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
def find_distance_between_triangles(triangle1, triangle2):
    distances1 = map(
            lambda x: (x[1][0], ([x[0]], x[1][1])),
            [(point, find_distance_to_triangle(triangle2, point)) for point in triangle1]
        )
    #print(distances1)
    distances2 = map(
            lambda x: (x[1][0], (x[1][1], [x[0]])),
            [(point, find_distance_to_triangle(triangle1, point)) for point in triangle2]
            )
    #print(distances2)
    distances3 = [
                x
                for line_segment1 in get_simplices_1(triangle1)
                for line_segment2 in get_simplices_1(triangle2)
                for x in find_distance_between_line_segments(line_segment1, line_segment2)
            ]
    #print(distances3)
    return select_by_some_value(distances1 + distances2 + distances3, min)
#method returns tuple - nearest distance and corresponding tuple of arrays of atoms from 1st and second tetrahedra
# returns [(<distance>, ([<points from surf1>], [<points from surf2>]))]
def find_distance(tetrahedra1, tetrahedra2):
    #print("in find distance")
    dist =  [
            x
            for triangle1 in get_simplices_1(tetrahedra1)
            for triangle2 in get_simplices_1(tetrahedra2)
            for x in find_distance_between_triangles(triangle1, triangle2)
        ]
    #print(dist)
    ret = select_by_some_value(dist, min)
    #print(ret)
    #print("after find distance")
    return ret
    #print("ok")

def get_tetrahedra(simplex, points) : return (map(lambda x: points[x], simplex))

#simple and super-slow: not actual hausdorrf, returns just pairs of nearest atoms now. and converts them to corresponding aa pairs
def hausdorff_distance(surface1, surface2):
    distances = [
        x
        for simplex1 in surface1[1]
        for simplex2 in surface2[1]
        for x in find_distance(get_tetrahedra(simplex1, surface1[0]), get_tetrahedra(simplex2, surface2[0]))
    ]
    #print("got dist")
    nearest_atoms = select_by_some_value(distances, min)
    def atoms_to_aa(atoms) :
        #print("atoms")
        #print(atoms[0].get_full_id())
        return list(set([a.get_parent().get_id()[1] for a in atoms]))
    #print("got nearest_atoms")
    #print(nearest_atoms)
    #res = map(lambda x: (
    #        x[0], (atoms_to_aa(x[1][0]), atoms_to_aa(x[1][1]))
    #        ), nearest_atoms)
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
    pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    print("got surfaces, next should find dist between them")
    dist = hausdorff_distance(surface1, surface2)
    print(dist)
