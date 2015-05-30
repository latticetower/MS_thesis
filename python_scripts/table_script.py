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
                'notFoundHotspots',
                'method' # method used to collect this data: string, either EHRA or other methods
                ])


import csv, optparse
from rcsb_loader import fetch_pdb
import numpy as np

def print_not_found(pdb, chain1, n1, n2, result3, chain2 = ''):
    output = open("results/ehra_{0}_{1}_{2}_not_found.txt".format(pdb, chain1, chain2), 'w')
    output.write("{0}: {1} - {2}\n ehra aa:\n".format(pdb, chain1, chain2))
    output.write("\t".join(map(str, n1))+"\nhotspots: \n")
    output.write("\t".join(map(str, n2))+"\nnot found:\n")
    output.write("\t".join(map(str, result3))+"\n")

    output.close()
if __name__ == "__main__":
    table_data = []
    unprocessed_data = dict()
    import logging, os
    logging.basicConfig(filename="logs/" + os.path.splitext(os.path.basename(__file__))[0] + ".log", level=logging.DEBUG)
    logging.debug("============\nCalled table script")
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
    file_names = dict()
    for pdb in names:
        file_names[pdb] = fetch_pdb(pdb, options.pdb_folder)
    from itertools import groupby
    tuples = sorted(tuples, key=lambda x: x.PDB_ID)
    for pdb, data in groupby(tuples, lambda x: x.PDB_ID):
        print("---------\n{}".format(pdb))
        from make_it_simple import read_pdb_info, main_func, read_pdb_chains_info
        from os import sep
        filename = file_names[pdb]
        l = list(data)
        logging.debug("-------------\nstart to process new pdb file")
        logging.debug(pdb)
        l = sorted(l, key=lambda x: x.Chain_ID)
        for chain, data in groupby(l, lambda x: x.Chain_ID):
            logging.debug("+++ chain N {}".format(chain))
            lchain = list(data)
            res = main_func(filename, chain, None)
            if res == None:
                logging.debug("main func not finished well")
                if not unprocessed_data.has_key(pdb): unprocessed_data[pdb] = dict()
                unprocessed_data[pdb][chain] = lchain
                continue
            n1 = np.array([int(y.get_id()[1]) for y in res[0] if y.get_id()[0]!='W'])
            #print(n1)
            n2 = np.array([int(x.PDB_res) for x in lchain])
            #print(n2)
            logging.debug(len(n1))
            logging.debug(len(n2))
            result = np.setdiff1d(n1, n2)
            #print(result)
            logging.debug("empty aa in selected region:{}".format(len(result)))
            result2 = np.intersect1d(n1, n2)
            logging.debug("Found hotspots: {}".format(len(result2)))
            result3 = np.setdiff1d(n2, n1)
            logging.debug("Not found hotspots: {}".format(len(result3)))
            logging.debug(result3)
            if len(result3) > 0:
                print_not_found(pdb, chain, n1, n2, result3)
            table_data.append(TableLine(pdb,
                chain, '-', res[1], len(n1),
                0, len(result2), len(n2), len(result3), 'ehra'))
    correct_output = open('ehra_results_2chains.txt', 'w')
    for tl in table_data:
        correct_output.write("\t".join(map(str, list(tl))) + "\n")
    correct_output.close()
    incorrect_output = open('morethan2chains.txt', 'w')
    for pdb in unprocessed_data.keys():
        incorrect_output.write(pdb + "\n")
    incorrect_output.close()
    table_data2 = []
    for pdb in unprocessed_data.keys():
        print(pdb)
        filename = file_names[pdb]
        chains = read_pdb_chains_info(filename)
        for chain1 in unprocessed_data[pdb].keys():
            lchain = unprocessed_data[pdb][chain1]
            total_hotspots = np.array([x.PDB_res for x in lchain])
            found_hotspots = np.asarray([])
            ehra_regions = dict()
            for chain2 in chains:
                if chain2 != chain1:
                    res = main_func(filename, chain1, chain2)
                    ehra_regions[chain2] = np.array([y.get_id()[1] for y in res[0] if y.get_id()[0]!='W'])
                    result3 = np.setdiff1d(ehra_regions[chain2], total_hotspots)
                    found_hotspots = np.union1d(found_hotspots, np.intersect1d(ehra_regions[chain2], total_hotspots))
                    table_data2.append(TableLine(pdb,
                        chain1, chain2, res[1], len(ehra_regions[chain2]),
                        0, 0, len(total_hotspots), len(result3), 'ehra'))
                    #print res
            table_data2.append(TableLine(pdb, chain1, '-',
                res[1], 0, 0,0, len(total_hotspots), len(found_hotspots), 'ehra total'))
            if len(found_hotspots) < len(total_hotspots):
                print_not_found(pdb, chain1, total_hotspots, np.unique(np.asarray([x for l in ehra_regions.values() for x in l])), found_hotspots, 'multi')
            #print data
        print chains
    correct_output = open('ehra_results_multichain.txt', 'w')
    for tl in table_data2:
        correct_output.write("\t".join(map(str, list(tl))) + "\n")
    correct_output.close()
    logging.debug("end of table script\n ---------------")
