
import numpy as np
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from ateams.structures import Lattice
from ateams.models import SwendsenWang
from ateams.stats import constant, critical
from ateams import Chain, Tape, _version
import sys

L = Lattice().fromCorners([8,8,8,8], dimension=2, field=3)

SW = SwendsenWang(L, temperatureFunction=constant(-0.585), LinBox=True)
N = 10
M = Chain(SW, steps=N)

for (spins, occupied) in M.progress():
	pass
