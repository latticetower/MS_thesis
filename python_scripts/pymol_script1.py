#!/usr/bin/env python

#the following script should load alascan data from asedb for given record and apply it to the structure
#it should also contain function for mesh creation with given cutoff value
#
#
from pymol import cmd, stored

def selectInterface(ch1, ch2, cutoff = 3.5):
    '''
DESCRIPTION
    some useful example scripts:
    #http://www.pymolwiki.org/index.php/FindObjectsNearby - that is useless, it returns object, not atoms or chains or smth
    http://www.pymolwiki.org/index.php/Get_Coordinates_I
    http://www.pymolwiki.org/index.php/DistancesRH
    http://www.pymolwiki.org/index.php/Find_pairs
    '''
    cmd.bg_color("white")
    cmd.select("standard_interface",
                "byres (({0} within {2} of {1}) or ({1} within {2} of {0}))".format(
                        ch1, ch2, float(cutoff)
                ))
    cmd.color("red", "standard_interface")
    #commented, because this can be useful
    #for ((m1, a1), (m2, a2)) in cmd.find_pairs(ch1, ch2, 1, 1, float(cutoff)):
    #    #print(1)
    #    atom = cmd.get_model(m1).atom[a1]
    #    #str = "%s/%s/%s/%s" % (atom.chain, atom.resn, atom.resi, atom.name)
    #    #print str
    #    cmd.select("residuedata", "or (chain {0} and resi {1})".format(atom.chain, atom.resi))
    #    #cmd.show('dot', "residuedata")
    #    cmd.color("red", "residuedata")
    #    #if atom.chain == "A":
    #    #cmd.select()
    #    #print m1
    #print stored.objs.keys() #return stored.objs.keys()
def build_selection_string(row, cursor, mut_obj_name, sysId, probable_chain = None):
    print ("in build_selection_string")
    sql_str = "SELECT AA, Chain, position, pos_suff, ddG FROM hotspot.mutation where sysID = {0}".format(sysId)
    cursor.execute(sql_str)
    d1 = cursor.fetchall()
    total_str = []
    from Bio.PDB.Polypeptide import one_to_three
    print row
    for (aa, chain, pos, pos_suff, ddG) in d1:
        print(ddG)
        if chain == '':
            if probable_chain == None:
                print("no chain specified, unimplemented")
                print(row)
                #todo: should show header info here from compound area
                return
            chain = probable_chain
        #here we select aminoacid for given chain and position
        if(abs(float(ddG)) >= 1.0):
            str = "(chain %s and resi %s and resn %s)" % (chain, pos, one_to_three(aa))
            print(str)
            total_str.append(str)
    return "byres (" + " or ".join(total_str) + ")"

# method shows and selects all aa with significant mutations from asedb if possible
# ignores and shows message if not possible
def selectMutants(mut_obj_name, sysId = None, probable_chain = None):
    #ignore if
    import MySQLdb
    import os
    db = MySQLdb.connect("localhost","root","","hotspot" )
    cursor = db.cursor()
    str = "SELECT * FROM hotspot.system where PDB like \"{0}\" ".format(mut_obj_name)
    if sysId != None:
        str += " and sysId = {0}".format(sysId)
    str += ";"
    cursor.execute(str)
    data = cursor.fetchall()
    if (len(data) > 1):
        print("Cannot decide what record to use, please select one and call selectMutants again with second parameter provided")
        for row in data:
            print row
        return
    if (len(data) < 1):
        print("Wrong params, exiting (no data obtained from database)")
        return
    # here we have exactly 1 record selected from hotspot.system
    str = (build_selection_string(data[0], cursor, mut_obj_name, sysId, probable_chain))
    print(str)
    cursor.close()
    r = cmd.select("mutants_{0}, ".format(sysId), str)
    print r
    cmd.color("blue", "mutants_{0}, ".format(sysId))
    #for row in data:
    #    trg_code = row[4]
    #    trg_file = pdb_dir + os.sep + trg_code.lower() +".pdb"
    #    process_pdb(row, trg_file, cursor, handle, handle2)
    #    #break
    db.close()

cmd.extend("selectInterface", selectInterface);
cmd.extend("selectMutants", selectMutants);

# how to call from pymol:
# 0. load pdb by hand (i.e. 3hhr.pdb)
# 1. run github/MS_thesis/python_scripts/pymol_script1.py
# 2. selectInterface A//, B//, 8
# 3. selectMutants 3hhr,1,A
# 4. in pymol "hide all", "show ((mutants_1) in (not standard_interface))" (with cutoff=8 shows nothing)
