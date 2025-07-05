
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import Bernoulli
import sys


try:
    DIM = int(sys.argv[-2])
    L = int(sys.argv[-1])
    sparse = False
except:
    L = 6
    DIM = 4

pyximport.install()

M = Bernoulli.construct(L, DIM)


TESTS = [L, DIM]
WIDTHS = [8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(30)


cProfile.runctx("Bernoulli.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/Bernoulli/{L}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

