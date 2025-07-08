
import warnings
warnings.simplefilter("ignore", UserWarning)

from Profile import profile
import IC
import sys


try:
    DIM = int(sys.argv[-3])
    FIELD = int(sys.argv[-1])
    L = int(sys.argv[-2])
except:
    FIELD = 3
    L = 10
    DIM = 4

M = IC.construct(L, DIM, FIELD)

TESTS = [L, DIM, f"Z/{FIELD}Z", bool(int(FIELD < 3)), M.model.faces]
WIDTHS = [8, 8, 8, 8, 8]
DESC = [str(thing).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)
DESC = ("      "+DESC).ljust(10)

profile(IC.chain, [M,DESC], f"./profiles/InvadedCluster/{L}.{DIM}.{FIELD}.{'PHAT.' if FIELD<3 else ''}txt")