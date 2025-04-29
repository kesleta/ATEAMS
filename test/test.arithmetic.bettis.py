
from ateams.arithmetic import Persistence, FINT
from ateams.structures import Lattice
from galois import GF
from math import comb as combinations
from random import choice
import numpy as np
import sys
from tqdm import tqdm

F = 3

# Number of trials.
SCALE = 3
DIMENSION = 4
LATTICES = [Lattice().fromCorners([SCALE]*d, field=F) for d in range(2, DIMENSION+1)]
SUBCOMPLEXES = 20

DESC = (f"          SPARSE, SERIAL").ljust(30)

with tqdm(total=len(LATTICES)*SUBCOMPLEXES, dynamic_ncols=True, desc=DESC) as bar:
	for lattice in LATTICES:
		for _ in range(SUBCOMPLEXES):
			P = Persistence(F, lattice.flattened)

			# First, check that we correctly compute the betti numbers for the full
			# torus.
			dimension = lattice.dimension
			ground = [combinations(dimension, d) for d in range(dimension+1)]

			print("######## BEGIN FULL COMPLEX COMPUTATION")
			subcomplex = np.array(range(len(lattice.flattened)))
			test = P.ComputeBettiNumbers(subcomplex)
			print("######## END FULL COMPLEX COMPUTATION")
			print()

			try:
				assert ground == test
				print(ground)
				print(test)
				print("#################################")
				print("#################################")
			except AssertionError:
				print(ground)
				print(test)
				print("BAD")
				print("#################################")
				print("#################################")
				exit(1)

			# Next, check what we get for a "random" subcomplex; we check that,
			# by deleting one cell of dimension i (for i > 1), we cannot realize
			# any cells of dimension higher than i.
			for d in range(1, dimension+1):
				delete = choice(range(*lattice.tranches[d]))
				subsubcomplex = np.concatenate([subcomplex[:delete], subcomplex[delete+1:]])
				ground = [combinations(dimension, d) for d in range(dimension+1)]

				print("######## BEGIN SUBCOMPLEX COMPUTATION")
				print("######## DIMENSION", d)
				print("######## DELETED", delete)
				print("######## SUBCOMPLEX")
				print("\t", subsubcomplex)
				print()
				test = P.ComputeBettiNumbers(subsubcomplex)

				try:
					assert ground != test
					print(ground)
					print(test)
				except AssertionError:
					print(ground)
					print(test)
					print("BAD")
					exit(1)
				
				print()
				print("######## END SUBCOMPLEX COMPUTATION")

			bar.update()

exit(0)