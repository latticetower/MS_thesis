#!/usr/bin/env python
import numpy as np

#method returns set of convex hull triangles from surface1, whose points lay within given cutoff near surface2
def find_interface_triangles(surface1, surface2, cutoff):
    triangles = surface1[1].convex_hull[np.any(surface2[2].query(
        surface1[1].points[surface1[1].convex_hull], 1, 0, 2, cutoff
        )[0] < cutoff, axis = 1)]
    return triangles


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
