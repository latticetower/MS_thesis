# #  _!/usr/bin/env python

import tkSimpleDialog, tkMessageBox, tkFileDialog
from pymol import cmd
import sys, urllib2, zlib

import pyximport; pyximport.install()
import logging, os, sys, inspect

cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
if cmd_folder not in sys.path:
    sys.path.insert(0, cmd_folder)
print(cmd_folder)

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0], "lib")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
print(cmd_subfolder)
import rcsb_loader

def __init__(self):
    self.menuBar.addmenuitem('Plugin', 'command',
                            'EHRA tools',
                            label = ' find energy hotspot regions',
                            command = lambda s=self : openFromFile(s))

from Tkinter import *
# plugin dialog class
class EHRADialog(tkSimpleDialog.Dialog):
    def body(self, master):
        self.parent = master
        Label(master, text="Click to select PDB file path:").grid(row=0)
        self.pdb_file_path = Button(master, text="select pdb file...", command=self.open_file).grid(row=0, column=1)#Entry(master)
        #self.pdb_file_path.bind("<Button-1>", self.open_file)
        #self.pdb_file_path.grid(row=0, column=1)
        #Button(master, text="RCSB", command=self.load_from_rcsb).grid(row=0, column=2)
        #optional button to enable file loading from remote rcsb ftp, currently turned off
        Label(master, text="Chain 1 (mutated):").grid(row=1)
        Label(master, text="Chain 2:").grid(row=2)
        self.ch1 = Entry(master)
        self.ch1.grid(row=1, column=1)
        self.ch1.insert(0, "L")
        # chain 2
        self.ch2 = Entry(master)
        self.ch2.insert(0, "H")
        self.ch2.grid(row=2, column=1)
        #dist cutoff
        Label(master, text="Distance cutoff:").grid(row=3)
        self.distance_cutoff = Entry(master)
        self.distance_cutoff.insert(0, "3.5")
        self.distance_cutoff.grid(row=3, column=1)
        #solvent radius
        Label(master, text="Solvent radius").grid(row=1)
        self.sas_radius = Entry(master)
        self.sas_radius.insert(0,"1.4")
        self.sas_radius.grid(row=4, column=1)
    def open_file(self, event = None):
        file_name = tkFileDialog.askopenfilename(parent=self.parent,
            title='EHRA-finder input PDB file selection')
        if file_name != None:
            print file_name
            self.pdb_file = file_name
        else:
            print 0
    #helper method. loads data from rcsb
    def load_from_rcsb(self, event=None):
        pdbCode = tkSimpleDialog.askstring('PDB ID to load',
            'Enter a 4-digit pdb code:',
            parent=self.parent)
        if pdbCode != None:
            saveDialogResult = tkFileDialog.asksaveasfilename(
                defaultextension = ".pdb",
                initialfile = pdbCode + ".pdb")
            from rcsb_loader import fetch_pdb
            fetch_pdb(pdbCode, os.path.dirname(saveDialogResult))
        #fetch_pdb(pdbCode)
    def apply(self):
        import string
        pdb_file = self.pdb_file
        if pdb_file == None:
            print "select pdb file!"
            return
        chain1 = self.ch1.get()
        chain2 = self.ch2.get() #string.atoi(
        sas_radius = string.atof(self.sas_radius.get())
        cutoff = string.atof(self.distance_cutoff.get())
        # load struct
        cmd.bg_color("white")
        cmd.load(pdb_file)
        cmd.h_add()
        cmd.remove("solvent")
        cmd.save(pdb_file)
        from ehra_finder import find_regions
        print "Processing, please wait..."
        print (chain1, chain2)
        res = find_regions(pdb_file, chain1, chain2, cutoff, sas_radius)
        str = "chain %s and (%s)" % ('L', " or ".join(["resi {0}".format(
                aa_info.get_id()[1]
                )
            for aa_info in res[0]]))
        # select region
        cmd.select("EHRA_region", "!(all)")
        cmd.select("EHRA_region", str)
        for chain in (chain1, chain2):
            cmd.select("chain_" + chain, "!(all)")
            cmd.select("chain_" + chain, chain + "//")
        cmd.hide("(all)")
        cmd.show("cartoon", "chain_" + chain1)
        cmd.show("lines", "chain_" + chain2)
        cmd.color("blue", "chain_" + chain1)
        cmd.color("green", "chain_" + chain1)
        cmd.color("red", "EHRA_region")
        cmd.save(os.path.dirname(pdb_file) + os.sep + chain1 + "_" + chain2 + ".wrl", "EHRA_region")
        print("done")
        # save to mesh file


# commands
#

def openFromFile(app):  EHRADialog(app.root)


if __name__ == "__main__":
    print "This script is a plugin for PyMOL and should be called from pymol gui"
