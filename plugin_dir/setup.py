#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize

setup(
  ext_modules = cythonize("lib/*.pyx"),
  package_dir = {'': 'bin'}
)
