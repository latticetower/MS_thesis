#!/usr/bin/env python
import pyximport; pyximport.install()
import logging, os, sys, inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

import rcsb_loader
if __name__ == "__main__":
    from ehra_finder import find_regions
    #logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    #logging.debug("============\nCalled simple script")
    res = find_regions(os.path.abspath('../test_data/2OSL.pdb'), 'L', 'H')
    #print res[0]
    str = "chain %s and (%s)" % ('L', " or ".join(["resi {0}".format(
            aa_info.get_id()[1]
            )
        for aa_info in res[0]]))
    print str
    print("end of processing")
