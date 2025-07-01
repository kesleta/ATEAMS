
import numpy as np
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from ateams.complex import Cubical
from ateams.model import Glauber
from ateams.stats import constant, critical
from ateams import Chain, Tape, _version
import sys

L = Cubical().fromCorners([3,3,3,3], field=3)

SW = Glauber(L, temperature=constant(-0.585))
N = 10
M = Chain(SW, steps=N)

for (spins, occupied) in M.progress():
	print(spins)
