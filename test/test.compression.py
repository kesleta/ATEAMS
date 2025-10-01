
from ateams import Chain, Recorder, Player
from ateams.statistics import critical
from ateams.complexes import Cubical
from ateams.models import Glauber, SwendsenWang
import numpy as np


C = Cubical().fromCorners([10,10,10,10])
model = SwendsenWang(C, field=7, dimension=2, temperature=critical(7))
M = Chain(model, steps=1000)

rec = Recorder()

actuals = []

with rec.record("data/out.lz", blocksize=23) as rec:
	for (spins, satisfied) in M.progress():
		rec.store((spins, satisfied))
		actuals.append((spins, satisfied))


play = Player()
saveds = []

with play.playback("data/out.lz", steps=1000) as play:
	i = 0
	for (spins, satisfied) in play.progress():
		actualspins, actualsatisfieds = actuals[i]
		assert (actualspins == spins).all()
		assert (actualsatisfieds == satisfied).all()
		i += 1