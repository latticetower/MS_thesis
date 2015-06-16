1. call from command line
'''
python setup.py build_ext --inplace
'''
2. after that, plugin dir (called 'ehra_tools') will contain compiled files.
gzip that folder.

3. after that, open PyMol, call Plugin -> Plugin Manager from main menu. In install dialog, select .zip file you previously created.

4. Now Plugin menu should contain new item - "find energy hotspots". This shows dialog window, specify PDB file, chain names and other parameters.
Than wait, and you'll see loaded data with specificity determining region.
