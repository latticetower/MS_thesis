#!/usr/bin/python

from Bio.PDB import *

cutoff = 5

pdb_dir = "pdb_folder"

import MySQLdb

max_pairs = 10
import os

c_hash = {
    'hGH' : 'human growth hormone',
    'RNase inhibitor': 'ribonuclease inhibitor',
    'Protein A': 'fragment b of protein a complex',
    'Factor VII': 'blood coagulation factor viia light chain',
    'hGHbp': 'human growth hormone receptor (hghbp)',
    'Im9': 'protein (colicin e9 immunity protein)',
    'BPTI': 'trypsin inhibitor'
}

def get_chain(compound, row1, row2):
    ones = []
    for x in compound:
        if len(compound[x]['chain'].split(',')) == 1:
            ones.append(x)
    #if len(ones) >= 1:
    #    return compound[ones[0]]['chain'].upper()
    for x in compound:
        if compound[x]['molecule'].upper() == row1.upper():
            #print row1.upper()
            #print compound[x]
            if len(compound[x]['chain'].strip(',')) > 1:
                print "cannot decide, there are more than 1 chain"
                return None
            return compound[x]['chain'].upper()
        if c_hash.has_key(row1):
            if compound[x]['molecule'].upper() == c_hash[row1].upper():
                return compound[x]['chain'].upper()
    print "cannot decide what compound to prefer"
    print compound
    print row1
    print row2
    return None
def get_chain_except(compound, ch):
    chains = set()
    for x in compound:
        if (ch not in map(lambda x: x.upper(), compound[x]['chain'].split(','))):
            for c in compound[x]['chain'].split(','):
                if c != ch:
                    chains.add(c)
    chains = map(lambda x: x.strip().upper(), chains)
    return chains

def  find_nearest_dist_to_other_chains(aa, c, structure, ddG, handle2):
    chains = get_chain_except(structure.header['compound'], c)
    dist = aa['CA'] - list(structure[0][chains[0]].get_atoms())[0]
    for chain in chains:
        #print chain
        for res in structure[0][chain]:
            for atom in res:
                for atom2 in aa:
                    dist = min(dist, atom - atom2)
    if dist > cutoff:
        return "{0} (ddG = {1}): {2}".format(dist, ddG, aa)
    return ""

def process_2chained(row, structure, filename, cursor, handle2):
    sql_str = "SELECT AA, Chain, position, pos_suff, ddG FROM hotspot.mutation where sysID = {0}".format(row[0])
    cursor.execute(sql_str)
    d1 = cursor.fetchall()
    for (aa, c, pos, pos_suff, ddG) in d1:
        #print row[0]
        if c == '':
            c = get_chain(structure.header['compound'], row[1], row[2])
            #print c
            if c == None: return
        if c in structure[0]:
            if len(structure[0][c]) >= pos:
                try:
                    str1 = find_nearest_dist_to_other_chains(structure[0][c][pos], c, structure, ddG, handle2)
                except KeyError:
                    print "got error"
                    str1 = ""
                if str1 != "":
                    handle2.write("\nat position {0}: \n{1}".format(pos, str1))
    print 1

def process_pdb(row, filename, cursor, handle, handle2):
    parser = PDBParser()
    structure = parser.get_structure(row[1] + '+' + row[2], filename)
    if len(structure.header['compound']) == 2:
        handle2.write('\n'+str(structure.header['compound']) + '\n')
        handle2.write("pdb = '{3}', sysId = {0}, mutated = '{1}', 2 = '{2}'\n".format(row[0], row[1], row[2], row[4]))
        process_2chained(row, structure, filename, cursor, handle2)
    else:
        handle.write("unprocessed:\n")
        handle.write( filename+'\n')
        handle.write( ', '.join(map(lambda x: str(x), row)))
        handle.write('\n')
        handle.write(str(structure.header['compound']))
        handle.write('\n------\n')
    print "next ->"


db = MySQLdb.connect("localhost","root","","hotspot" )
cursor = db.cursor()

cursor.execute("SELECT * FROM hotspot.system where PDB!=0;")

data = cursor.fetchall()




handle = open("unprocessed.txt",'w')
handle2 = open("mutated.txt", 'w')
for row in data:
    trg_code = row[4]
    trg_file = pdb_dir + os.sep + trg_code.lower() +".pdb"
    process_pdb(row, trg_file, cursor, handle, handle2)
    #break
handle2.close()
handle.close()
db.close()
