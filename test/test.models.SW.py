
import numpy as np
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from ateams.structures import Lattice
from ateams.models import SwendsenWang, CPSwendsenWang
from ateams.stats import constant, critical
from ateams import Chain, Tape, _version
import sys

L = Lattice().fromCorners([4,4,4,4], field=3)

SW = CPSwendsenWang(L, temperatureFunction=constant(-0.585))
N = 25
M = Chain(SW, steps=N)

for (spins, occupied) in M.progress():
	print(spins)
	print(occupied)
	print()
