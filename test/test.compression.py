
from ateams import Chain, Recorder, Player
from ateams.statistics import critical
from ateams.complexes import Cubical
from ateams.models import Glauber, SwendsenWang
import numpy as np


C = Cubical().fromCorners([10,10,10,10])
model = SwendsenWang(C, field=5, dimension=2, temperature=critical(5))
M = Chain(model, steps=1000)

actuals = []

with Recorder().record("data/out.lz", blocksize=23) as rec:
	for (spins, satisfied) in M.progress():
		rec.store((spins, satisfied))


saveds = []

with Player().playback("data/out.lz", steps=1000) as play:
	for (spins, satisfied) in play.progress():
		pass