
import numpy as np
import dateutil.relativedelta, time, datetime, json, pathlib, sys, platform
from ateams.complexes import Cubical
from ateams.models import Glauber
from ateams.statistics import constant, critical
from ateams import Chain, Tape, _version
import sys

L = Cubical().fromCorners([10,10,10,10], field=3)

SW = Glauber(L, temperature=constant(-0.585))
N = 10
M = Chain(SW, steps=N)

for (spins, occupied) in M.progress():
	pass
