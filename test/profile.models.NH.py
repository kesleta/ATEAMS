
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import NH
import sys


try:
    field = int(sys.argv[-1])
    L = int(sys.argv[-2])
except:
    field = 3
    L = 10

pyximport.install()

M = NH.construct(L, field)

TESTS = [L, f"Z/{field}Z"]
WIDTHS = [8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)

cProfile.runctx("NH.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/Nienhuis/{L}.{field}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

