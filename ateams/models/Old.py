 
import numpy as np

from ..arithmetic import MatrixReduction, KernelSample, FINT, evaluateCochain
from ..stats import constant
from .Model import Model


class Nienhuis(Model):
	name = "Nienhuis"
	
	def __init__(self, 
			L, q1, q2, temperatureFunction=constant(-0.6), initial=None,
			sparse=True, parallel=False, minBlockSize=32, maxBlockSize=64, cores=4
		):
		"""
		Initializes the self-dual Nienhuis model.

		Args:
			L: The `Lattice` object on which we'll be running experiments.
			q1 (float): Edge coupling parameter.
			q2 (float): Face coupling parameter.
			temperatureFunction (Callable): A temperature schedule function which
				takes a single positive integer argument `t`, and returns the
				scheduled temperature at time `t`.
			initial (galois.FieldArray): A vector of spin assignments to components.
			sparse (boolean): Should matrices be formatted sparsely? (Uses C/C++).
			parallel (boolean): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.
		"""

		self.lattice = L
		self.temperatureFunction = temperatureFunction #DOESNT DO ANYTHING
		self.faceCells = len(self.lattice.boundary[self.lattice.dimension-1])
		self.cubeCells = len(self.lattice.boundary[self.lattice.dimension])
		self.q1=q1
		self.q2=q2
		self.p1=q1/(1+q1)
		self.p2=q2/(1+q2)

		self.sparse = sparse
		self.parallel = parallel
		self.minBlockSize = minBlockSize
		self.maxBlockSize = maxBlockSize
		self.cores = cores

		# SW defaults.
		self.spins = initial if initial else self.initial()
		
		self.coboundary = self.lattice.matrices.coboundary.astype(FINT)
		self.Reducer = MatrixReduction(self.lattice.field.characteristic, parallel, minBlockSize, maxBlockSize, cores)
		self.SampleFromKernel = KernelSample


	def initial(self):
		"""
		Computes an initial state for the model's Lattice.

		Returns:
			A Galois `Array` representing a vector of spin assignments.
		"""
		return self.lattice.field.Random(self.faceCells)
	

	def proposal(self, time):
		"""
		Proposal scheme for generalized Swendsen-Wang evolution on the Potts model.

		Args:
			time (int): Step in the chain.

		Returns:
			A Galois `Array` representing a vector of spin assignments.
		"""
		# Compute the probability of choosing any individual cube in the complex.
		p1=self.p1
		p2=self.p2
		
		# Choose cubes to include; in effect, this just does a boatload of indexing.
		uniformsFace = np.random.uniform(size=self.cubeCells) #randomly assigns a probability to each face
		uniformEdge = np.random.uniform(size=self.faceCells) #randomly assigns a probability to each edge
		includeFace = (uniformsFace < p1).nonzero()[0] #gives you indices where probability is under threshold
		includeEdge =(uniformEdge < p2).nonzero()[0] #gives you indices where probability is under threshold

		boundary = self.lattice.boundary[self.lattice.dimension][includeFace] #list of edges for faces under threshold
		boundaryValues = evaluateCochain(boundary, self.spins) #list of face sums for open faces
		percedFaces = (boundaryValues == 0).nonzero()[0] #indices where is 0 (that is also selected by perc)
		
		percedEdges = np.intersect1d((self.spins == 0).nonzero()[0], includeEdge)
		includedEdges=np.setdiff1d(range(0,len(self.lattice.boundary[self.lattice.dimension-1])),percedEdges)
		

		#satisfied = np.zeros(len(self.lattice.boundary[self.lattice.dimension])).astype(int) #all 0's for faces
		#satisfied[zeroFaces] = 1 #all open faces which evaluate to 0 get label 1

		# Uniformly randomly sample a cocycle on the sublattice admitted by the
		# chosen edges; reconstruct the labeling on the entire lattice by
		# subbing in the values of c which differ from existing ones.
		subbasis = self.coboundary.take(percedFaces, axis=0).take(includedEdges, axis=1)
		changedSpins = self.SampleFromKernel(self.Reducer, subbasis)
		# changedSpins=sampleFromKernel(self.lattice.matrices.coboundary, self.lattice.field, relativeCells=percedFaces, relativeFaces=includedEdges) #based on this face set, find a cocycle, return cocycle and labelling

		newSpins=self.spins.copy()
		newSpins[includedEdges] = changedSpins
		newSpins[percedEdges] = 0

		edgeEnergy=(len(newSpins)-sum((newSpins==0)))/len(newSpins)

		fullBoundary = self.lattice.boundary[self.lattice.dimension]
		fullBoundaryValues=evaluateCochain(fullBoundary, self.spins)
		faceEnergy=(len(fullBoundaryValues)-sum((fullBoundaryValues==0)))/len(fullBoundaryValues)

		#energy=edgeEnergy+faceEnergy
		#print(edgeEnergy,faceEnergy)
		
		return newSpins, (edgeEnergy, faceEnergy)

	def assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (galois.FieldArray): Cocycle on the sublattice.
		
		Returns:
			None.
		"""
		self.spins = cocycle