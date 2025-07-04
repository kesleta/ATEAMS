
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import Glauber
import sys


try:
    L = int(sys.argv[-1])
    sparse = False
except:
    L = 10

pyximport.install()

M = Glauber.construct(L)

TESTS = [L]
WIDTHS = [8]
DESC = [str(int(thing)).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)


cProfile.runctx("Glauber.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/Glauber/{L}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

