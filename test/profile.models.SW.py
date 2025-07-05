
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import SW
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

M = SW.construct(L, dim, field)

TESTS = [L, dim, f"Z/{field}Z"]
WIDTHS = [8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)

FAIL = False
EXCEPT = None

try:
	cProfile.runctx("SW.chain(M, DESC)", globals(), locals(), "Profile.prof")
except Exception as e:
	EXCEPT = e
	FAIL = True

fname = f"./profiles/SwendsenWang/{dim}.{L}.{field}.txt"

with open(fname, 'w') as stream:
	if FAIL:
		stream.write(str(EXCEPT))
		exit(1)
	else:
		s = pstats.Stats('Profile.prof', stream=stream)
		s.strip_dirs().sort_stats("time").print_stats()

