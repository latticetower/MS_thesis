#!/usr/bin/env python

from triangulator import *
#import ala_scan

if __name__ == "__main__":
    pair_of_chains = read_pdb_info('../test_data/2OSL.pdb')
    surface1 = process_chain(pair_of_chains[0])
    surface2 = process_chain(pair_of_chains[1])
    print("got surfaces, next should find dist between them")
    dist = hausdorff_distance(surface1, surface2)
    print(dist)
