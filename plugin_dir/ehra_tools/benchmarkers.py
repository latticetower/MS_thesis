#!/usr/bin/env python

from Bio.PDB import *
import numpy as np
import timeit
import logging
import os

########################################
# decorator methods (for benchmarking) #
########################################
function_timings = {}

def timed(fn):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        result = fn(*args, **kwargs)
        elapsed = timeit.default_timer() - start_time
        logging.debug("Function \"" + fn.__name__ + "\" timed: " + str(elapsed))
        return result
    return wrapper

def accumulated(fn):
    def wrapper(*args, **kwargs):
        start_time = timeit.default_timer()
        result = fn(*args, **kwargs)
        elapsed = timeit.default_timer() - start_time
        if (function_timings.has_key(fn.__name__)):
            function_timings[fn.__name__][0] += 1
            function_timings[fn.__name__][1] += elapsed
        else:
            function_timings[fn.__name__] = [1, elapsed]
        #logging.debug("Function \"" + fn.__name__ + "\" timed: " + str(elapsed))
        return result
    return wrapper

import atexit


@atexit.register
def print_accumulated():
    logging.debug("Start of accumulated timings")
    for func_name in function_timings.keys():
        [total_calls, total_time] = function_timings[func_name]
        logging.debug(func_name + " : " + str(total_calls) + " calls, elapsed: " + str(total_time) +
            " average: " + str(total_time/total_calls))
    logging.debug("End of accumulated timings")
