#!/usr/bin/env python

from pymol import cmd

__all__ = [
    "alpha_shapes",
    "data_readers",
    "benchmarkers",
    "ehra_finder",
    "ehra_interface",
    "ehra_loops",
    "ehra_pockets",
    "ehra_tools",
    "rcsb_loader",
    "utils"
]

def make_global():
    '''
    "psico" might be installed as submodule of something else (PyMOL plugin).
    Invoke this function if you want to do "from psico import ...".
    '''
    import sys
    if sys.modules.get('psico') != sys.modules[__name__]:
        sys.modules['psico'] = sys.modules[__name__]

def init(pymolapi = 0):
    import pymol
    from pymol import cmd
    # init all submodules
    ehra_tools = __import__(__name__, fromlist=__all__)
    # pymol namespace
    if not hasattr(pymol, 'ehra_tools'):
        pymol.ehra_tools = ehra_tools
    if pymolapi:
        init_cmd(pymolapi == 2)

# PyMOL Plugin hook
def __init_plugin__(self=None):
    init(0)
    make_global()
    from ehra_tools import openFromFile
    self.menuBar.addmenuitem('Plugin', 'command',
                            'EHRA tools',
                            label = ' find energy hotspot regions',
                            command = lambda s=self : openFromFile(s))

def init_cmd(force=0):
    '''
    Adds all functions to PyMOL API (pymol.cmd).
    If force is True, overwrite existing names.
    '''
    from pymol import cmd
    for name, value in cmd.keyword.iteritems():
        function = value[0]
        if function.__module__.startswith(__name__):
            if force or not hasattr(cmd, function.__name__):
                setattr(cmd, function.__name__, function)
