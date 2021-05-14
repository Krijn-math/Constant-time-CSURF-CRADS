from progress.bar import Bar
from functools import reduce
import sys
import random
from math import ceil, floor, log, sqrt, pi
import numpy
import statistics
import argparse

# --------------------------------------------------------------------------------------------------------------------------------
def getinputs(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-a", "--algorithm", help="csidh or csurf algorithm", required=True)
    parser.add_argument("-p", "--prime", help="prime number configuration should be stored in pSUFFIX (sop folder is taken as default).", required=True)
    parser.add_argument("-m", "--multieval", help="unscaled or scaled multi-evaluation procedure", required=True)
    parser.add_argument("-s", "--style", help="style to be used: wd1 (with dummy operations and a single torsion point), wd2 (with dummy operations and a two torsion point), or df (dummy-free approach).")
    parser.add_argument("-b", "--benchmark", type=int, help="number of experiments to be used in the benchmark.", default=128)
    parser.add_argument("-e", "--exponent", type=int, help="For determining the number k of small odd primes to be used. The keyspace size is either (2e + 1)^k [wd2 style] or (e + 1)^n [wd1 and df styles].", default=1)
    parser.add_argument("-u", "--units", dest='units', action='store_true', help="Used to precompute the running time of each small odd prime degree isogeny construction and evaluation with velusqrt formulae (required only at first run).")
    parser.add_argument("-v", "--verbose", dest='verbose', action='store_true', help="Verbose mode.")
    parser.add_argument("-r", "--radical", dest='radicals', action='store_true', help="Radical isogenies mode (not only degree-2 isogenies on the surface.")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    options = parser.parse_args(args)
    return options

# Inputs
setting = getinputs(sys.argv[1:])
