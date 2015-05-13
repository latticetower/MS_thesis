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
                'totalHotspots'
                ])
