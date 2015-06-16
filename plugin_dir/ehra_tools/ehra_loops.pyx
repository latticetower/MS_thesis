#!/usr/bin/env python
from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

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
