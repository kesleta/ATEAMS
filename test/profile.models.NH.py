
import pstats, cProfile
import pyximport
import NH
import sys


try:
    L = int(sys.argv[-7])
    sparse = bool(int(sys.argv[-6]))
    parallel = bool(int(sys.argv[-5]))
    minBlockSize = int(sys.argv[-4])
    maxBlockSize = int(sys.argv[-3])
    cores = int(sys.argv[-2])
    slurm = bool(int(sys.argv[-1]))
except:
    L = 3
    sparse = False
    parallel = False
    minBlockSize = 32
    maxBlockSize = 64
    cores = 2
    slurm = False

pyximport.install()

M = NH.construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores)

DESC = [str(int(thing)).ljust(10) for thing in [L, sparse, parallel, minBlockSize, maxBlockSize, cores]]
DESC = " ".join(DESC)

cProfile.runctx("NH.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/Nienhuis/{L}{'.sparse' if sparse else ''}{'.parallel' if parallel else ''}{'.slurm' if slurm else ''}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

