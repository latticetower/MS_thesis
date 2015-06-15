#!/usr/bin/python

def fetch_pdb(pdbCode, pdb_dir):
    import urllib, gzip, os, string
    from os import sep
    outFile = pdb_dir + sep + pdbCode + ".pdb"
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
    return outFile

if __name__ == "__main__":
    pdb_dir = "pdb_folder"
    max_pairs = 10
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
    for row in data:
        trg_code = row[0]
        fetch_pdb(trg_code, pdb_dir)
