
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import SW
import sys


try:
    LinBox = int(sys.argv[-5])
    L = int(sys.argv[-4])
    parallel = bool(int(sys.argv[-3]))
    cores = int(sys.argv[-2])
    slurm = bool(int(sys.argv[-1]))
    sparse = False
except:
    LinBox = True
    L = 10
    parallel = False
    cores = True
    slurm = False
    sparse = False

pyximport.install()

M = SW.construct(L, parallel, cores, LinBox)

TESTS = [L, parallel, cores, LinBox]
WIDTHS = [8, 8, 8, 8]
DESC = [str(int(thing)).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)


cProfile.runctx("SW.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/SwendsenWang/{L}{'.sparse' if sparse else ''}{'.parallel' if parallel else ''}{'.slurm' if slurm else ''}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

