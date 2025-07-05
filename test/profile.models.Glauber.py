
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import Glauber
import sys


try:
    field = int(sys.argv[-3])
    DIM = int(sys.argv[-2])
    L = int(sys.argv[-1])
except Exception as e:
    L = 10
    DIM = 4
    field = 3

pyximport.install()

M = Glauber.construct(L, field, DIM)

TESTS = [L, DIM, f"Z/{field}Z"]
WIDTHS = [8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(20)

pr = cProfile.Profile()
pr.enable()
code = Glauber.chain(M, DESC)
pr.disable()

fname = f"./profiles/Glauber/{L}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

exit(code)

