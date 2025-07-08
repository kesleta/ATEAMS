
import warnings
warnings.simplefilter("ignore", UserWarning)

from Profile import profile
import Glauber
import sys


try:
    FIELD = int(sys.argv[-3])
    DIM = int(sys.argv[-2])
    L = int(sys.argv[-1])
except Exception as e:
    L = 10
    DIM = 4
    FIELD = 3


M = Glauber.construct(L, FIELD, DIM)

TESTS = [L, DIM, f"Z/{FIELD}Z", M.model.faces]
WIDTHS = [8, 8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(20)

profile(Glauber.chain, [M,DESC], f"./profiles/Glauber/{L}.{DIM}.{FIELD}.txt")
