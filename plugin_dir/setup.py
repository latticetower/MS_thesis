#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize

setup(
  ext_modules = cythonize("ehra_tools/*.pyx"),
  #package_dir = {'': 'bin'}
  name='ehra_finder',
  description='Energy Hotspot Regions selection tool',
  license='BSD-3-Clause',
  requires=['pymol'],
)
