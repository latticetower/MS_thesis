#!/usr/bin/env python

#the following script should select regions with hotspot residues and show them
#optionally to save them
#
# TODO: process correctly when chain doesn't present in structure
#
#
from pymol import cmd, stored
import sys, os, inspect

sys.path.append(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
import make_it_simple

# method selects aminoacids in energy hotspot residues area (EHRA)
# (some nice abbreviation wouldn't hurt anyone =))
def selectEHRA(ch1, ch2, filename, cutoff = 3.5):
    '''
DESCRIPTION
    Brief description what this function does goes here
    '''
    #
    # Your code goes here
    #
    #print(dir(make_it_simple))
    cutoff = float(cutoff)
    cmd.h_add(ch1)
    cmd.h_add(ch2)
    atoms1 = cmd.get_model(ch1).atom
    atoms2 = cmd.get_model(ch2).atom
    surface1 = process_atom_chain(list(atoms1))
    surface2 = process_atom_chain(list(atoms2))
    chains_ss_info = read_dssp_info(filename, ch1[0], ch2[0]) #this is temporary
    triangles = find_interface_triangles(surface1, surface2, cutoff)
    str = " or ".join(["(chain {0} and resi {1} and resn {2} and name {3})".format(
            aa_info.chain, aa_info.resi, aa_info.resn, aa_info.name)
        for aa_info in to_aa2(triangles, surface1)])
    str = " (%s)" % str
    cmd.select("EHRA_interface", "!(all)")
    cmd.select("EHRA_interface", str)
    from chempy import cpv
    check_edge = lambda x1, x2 : x1.vdw + x2.vdw < cpv.distance(x1.coord, x2.coord)
    check_edge2 = lambda x1, x2 : True
    def check_edge3(x1, x2):
        #print("{0} <> {1} ({2}, {3}) ".format(cpv.distance(x1.coord, x2.coord), (x1.vdw + x2.vdw + 1.4*2),x1.coord))
        dist = cpv.distance(x1.coord, x2.coord) - (x1.vdw + x2.vdw + 1.4*2)
        #print(dist)
        return dist > 0
    G = DTGraph(surface1)
    #nodes = G.find_pockets(triangles, check_edge3)
    #print(nodes)
    def distance_func(x1, x2):
        from chempy import cpv
        return cpv.distance(x1.coord, x2.coord)
    nodes = np.setdiff1d(G.find_pockets(triangles, check_edge3, distance_func), triangles)
    #cmd.select("EHRA_interface", )
    #print(to_aa2(triangles, surface1))
    str = "(%s)" % " or ".join(["(chain {0} and resi {1} and resn {2} and name {3})".format(
            aa_info.chain,
            aa_info.resi,
            aa_info.resn,
            aa_info.name
            )
        for aa_info in to_aa2(nodes, surface1)])
    cmd.select("EHRA_triang", str)
    print(chains_ss_info.keys())
    aa_with_coils = extend_to_coils(triangles, chains_ss_info[ch1[0]], surface1, lambda x: x.resi)
    print aa_with_coils
    str = "(%s)" % " or ".join(["(chain {0} and resi {1})".format(
            ch1[0],
            aa_info
            )
        for aa_info in aa_with_coils])
    print str
    cmd.select("EHRA_coils", str)
    #print(surface1)
    print "Hello, PyMOLers"
    print "You passed in %s and %s" % (ch1, ch2)
    print "I will return them to you in a list.  Here you go."
    return (ch1, ch2)

cmd.extend("selectEHRA", selectEHRA)
