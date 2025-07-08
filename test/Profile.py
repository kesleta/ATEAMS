
import pstats, cProfile
import sys
import pyximport
from pathlib import Path

def profile(func, args, fname):
	pyximport.install()
	pr = cProfile.Profile()
	pr.enable()
	code = func(*args)
	pr.disable()

	fname = Path(fname)
	fname.parent.mkdir(exist_ok=True, parents=True)

	with open(fname, 'w') as stream:
		s = pstats.Stats(pr, stream=stream)
		s.strip_dirs().sort_stats("time").print_stats()

	sys.exit(code)
