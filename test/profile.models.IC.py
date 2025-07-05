
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import IC
import sys


try:
    dim = int(sys.argv[-3])
    field = int(sys.argv[-1])
    L = int(sys.argv[-2])
except:
    field = 3
    L = 10
    dim = 4

pyximport.install()

M = IC.construct(L, dim, field)

TESTS = [L, dim, f"Z/{field}Z", bool(int(field < 3))]
WIDTHS = [8, 8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)


pr = cProfile.Profile()
pr.enable()
code = IC.chain(M, DESC)
pr.disable()

fname = f"./profiles/InvadedCluster/{L}.{field}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()


exit(code)
