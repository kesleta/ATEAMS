
import numpy as np
from ateams.complex import Cubical
from ateams.stats import constant, critical
from ateams import Chain, Tape, _version
import sys

L = Cubical().fromCorners([3,3], field=3)
L.toFile("lattice.json")

M = Cubical()
M.fromFile("lattice.json")
