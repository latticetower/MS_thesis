#!/usr/bin/python

pdb_dir = "pdb_folder"

max_pairs = 10

def fetch_pdb(pdbCode,outFile):
    import urllib
    import gzip
    import os
    import string

    remoteCode = string.upper(pdbCode)
    if not os.path.exists(pdb_dir):
        os.mkdir(pdb_dir)
    if not os.path.exists(outFile):
        try:
            filename = urllib.urlretrieve(
                'http://www.rcsb.org/pdb/cgi/export.cgi/' +
                remoteCode + '.pdb.gz?format=PDB&pdbId=' +
                remoteCode + '&compression=gz')[0]
        except:
            print "warning: %s not found.\n"%pdbCode
        else:
            if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
                try:
                    abort = 0
                    open(outFile, 'w').write(gzip.open(filename).read())
                    print "fetched: %s"%(pdbCode)
                except IOError:
                    abort = 1
                if abort:
                    os.remove(outFile)
            else:
                print "warning: %s not valid.\n"%pdbCode
            os.remove(filename)


import MySQLdb
db = MySQLdb.connect("localhost","root","","hotspot" )
# prepare a cursor object using cursor() method
cursor = db.cursor()

# execute SQL query using execute() method.
cursor.execute("SELECT distinct PDB FROM hotspot.system where PDB!=0;")

# Fetch a single row using fetchone() method.
data = cursor.fetchall()


# disconnect from server
db.close()
import os

for row in data:
    trg_code = row[0]
    trg_file = pdb_dir+os.sep+trg_code+".pdb"
    fetch_pdb(trg_code,trg_file)
