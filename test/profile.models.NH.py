
import warnings
warnings.simplefilter("ignore", UserWarning)

from Profile import profile
import NH
import sys


try:
    FIELD = int(sys.argv[-1])
    L = int(sys.argv[-2])
except:
    FIELD = 3
    L = 10

M = NH.construct(L, FIELD)

TESTS = [L, 2, f"Z/{FIELD}Z", M.model.faces]
WIDTHS = [8, 8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)

profile(NH.chain, [M,DESC], f"./profiles/Nienhuis/{L}.2.{FIELD}.txt")

