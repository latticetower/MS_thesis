#!/usr/bin/env python

# main script for final result tables generation
#


# method should return list\tuple of values:
# 1. pdb file identifier (string)
# 2. mutated chain
# 3. second chain (secondChain)
# 4. mutated chain length (mutatedLength)
# 5. number of aminoacids in hotspot regions in mutated chain (ehraLength)
# 6. number of masked aminoacids
# 7. number of hotspots found
# 8. total number of hotspots
# 9. number of hotspots found by interface scanning (with given threshold)
# 10. threshold
from collections import namedtuple
TableLine = namedtuple('TableLine', [
                'pdb',
                'mutatedChain',
                'secondChain',
                'mutatedLength',
                'ehraLength', # number of aminoacids in hotspot regions in mutated chain
                'maskedLength', # number of masked aminoacids
                'foundHotspots',
                'totalHotspots',
                'method' # method used to collect this data: string, either EHRA or other methods
                ])


import csv, optparse
from rcsb_loader import fetch_pdb

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('--csv', dest = 'csv', default = "",
        help = 'the csv file containing mutations information')
    parser.add_option('--pdb_folder', dest="pdb_folder", default="temp_pdb",
        help = "folder for saving pdb files from rcsb server")
    (options,args) = parser.parse_args()
    with open(options.csv) as infile:
        reader = csv.DictReader(infile)
        reader.fieldnames[:] = map(lambda x: x.translate(None, '#:(),'), reader.fieldnames)
        print(reader.fieldnames)
        Data = namedtuple('Data', reader.fieldnames)
        tuples = [Data(**row) for row in reader]
    names = set(map(lambda data:data.PDB_ID, tuples))
    for pdb in names:
        fetch_pdb(pdb, options.pdb_folder)
