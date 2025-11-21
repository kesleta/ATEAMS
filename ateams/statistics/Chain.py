
from . import always
from ..common import NumericalInstabilityWarning

import lz4.frame
import numpy as np


class Chain:
	"""
	Simulates a Markov chain on the given Model. For example, to generate 1000 samples
	from 2-dimensional Potts lattice gauge theory with coefficients in \(\mathbb Z/5\mathbb Z\)
	and store the spins on each 1-cell and the occupied 2-cells, use something like

	```python
	from ateams import Chain, Recorder
	from ateams.statistics import critical
	from ateams.complexes import Cubical
	from ateams.models import SwendsenWang

	# Construct the cubical complex; instantiate model, chain.
	C = Cubical().fromCorners([10]*4)
	SW = SwendsenWang(C, dimension=2, field=5, temperature=critical(5))
	M = Chain(SW, steps=1000)

	# Run the chain, storing the output.
	with Recorder().record("out.lz") as rec:
		for (spins, occupied) in M.progress():
			# You can store as many outputs as you'd like --- even ones
			# that aren't spit out by the model --- but they must be passed
			# to the `.store()` method as a tuple.
			rec.store((spins, occupied))
	```

	The `Recorder` class stores the differences between consecutive iterations using
	the [LZ4 compression schema](github.com/lz4/lz4), compressing the diffs in
	chunks of (default) size 100. You can replay the data with the `Player`
	class:
	
	```python
	from ateams import Player

	with Player().playback("out.lz", steps=1000) as play:
		for (spins, occupied) in play.progress():
			<do whatever>
	```

	Running the sampler and recording the data takes ~26 seconds (~38 it/s) on
	an M2 MacBook Air; replaying the data takes ~0 seconds (~24,211 it/s). The
	size of `out.lz` is ~1.8MB, so storing each cell's data requires \(1/11\)th
	of a byte per iteration (amortized).
	"""
	def __init__(self, model, accept=always(), statistics={}, steps=10000):
		"""
		Args:
			model (Model): A Model (e.g. `ateams.models.SwendsenWang`)
				generating the Markov chain.
			accept (Callable): A function that consumes the complex, model, and
				state to determine whether we're going to a good place.
			statistics (dict): A mapping of names to functions which take the complex,
				model, and state as argument. The Chain keeps track of these at
				each iteration and stores whatever output is given.
			steps (int): The number of iterations in the chain.
		"""
		self.model = model
		self.steps = steps
		self.accept = accept
		self._exitcode = 0
		self._warnings = 0

		# Store stats and things.
		self.functions = statistics
		self.statistics = { name: [] for name in self.functions.keys() }


	def __iter__(self):
		"""
		Initializes the Chain object as a generator.
		"""
		self.step = 0
		self.state = tuple([self.model.spins] + []*(self.model._returns-1))
		return self
	

	def __next__(self):
		"""
		Performs the computations specified by the proposal and acceptance schemes.
		"""
		# While we haven't reached the max number of steps, propose a new plan,
		# check whether it's acceptable/valid, and continue.
		while self.step < self.steps:
			# Propose the next state and check whether we want to accept it as
			# the next state or not; assign whichever state is chosen to the
			# Model.
			try: proposed = self.model._proposal(self.step)
			except NumericalInstabilityWarning:
				self._exitcode = 2
				self._warnings += 1
				proposed = self.state

			self.state = proposed if self.accept(self.state, proposed, self.step) else self.state
			self.model._assign(self.state[0])

			# Compute statistics.
			for name, function in self.functions.items():
				self.statistics[name].append(function(self.model, self.state))
			
			# Iterate.
			self.step += 1
			
			return self.state
		
		raise StopIteration
	
	
	def progress(self, dynamic_ncols=True, desc=""):
		"""
		Progress bar.

		Returns:
			`tqdm` iterable.
		"""
		from tqdm.auto import tqdm
		import sys

		# stream = sys.stdout if self.model._DEBUG else sys.stderr
		stream = sys.stderr
		return tqdm(self, file=stream, total=self.steps, dynamic_ncols=dynamic_ncols, desc=desc)


class Recorder:
	"""
	Safely "records" states from a Chain by tracking changes between successive
	states, and writing these changes --- alongside the initial state --- to
	file. WARNING: all input is flattened when recording and is NOT reshaped when
	playing back.
	"""
	def __init__(self): pass

	def record(self, fp: str, blocksize=100):
		"""
		Called to configure the recording apparatus for the Chain, and *should*
		be used as a context manager (i.e. with the `with` statement).

		Args:
			fp (str): Filename.
			block (int=100): Number of states to bundle together before compressing. 

		Returns:
			This `Recorder` object.

		
		"""
		# Configure filepath, writer.
		self._fp = fp
		self._writer = lz4.frame.LZ4FrameFile(f"{self._fp}", mode="wb", compression_level=9)
		
		# Save the chain for access during iteration; set the "previous" state to
		# be all -1s, since the first state yielded by iterating over the Chain
		# is the initial state, and we need to record the whole thing; fix a list
		# of possible integer states for later.
		self.previous = None
		self._previous = False

		self._blocksize = blocksize
		self._block = []
		self._stored = 0

		# Compression markers.
		self._intrablockbreak = '<<#>>'
		self._interblockbreak = '<<<##>>>'
		self._interlistbreak = '<>'
		
		# Enter the context manager.
		return self.__enter__()


	def store(self, state) -> None:
		"""
		Stores a state yielded by iteration over a `Chain`.

		Args:
			state (tuple): Tuple of NumPy arrays to be written to file.
		"""
		# Check whether we've received any arguments before; if not, initialize
		# a list to keep track of "previous" arguments. We expect the arguments
		# here to be numpy arrays.

		# Check whether we have any previously-stored data.
		if not self._previous:
			self.previous = tuple(-np.ones(shape=s.flatten().shape, dtype=s.dtype) for s in state)
			self._previous = True

		subblock = []
		
		# Get diffs; store.
		for i in range(len(state)):
			# Compute the diff and store.
			curr = state[i].flatten()
			diff = np.flatnonzero(~np.equal(curr, self.previous[i]))
			self.previous[i][diff] = curr[diff][::]

			# Add the encoded dictionary of diffs to the subblock.
			subblock.append(f"{' '.join(diff.astype(str))}{self._interlistbreak}{' '.join(curr[diff].astype(str))}")

		# Increment the number of iterations we've stored; if we've reached the
		# limit, compress, write to file, and begin again.
		self._block.append(subblock)
		self._stored += 1

		if self._stored == self._blocksize:
			self._writer.write(
				(self._interblockbreak.join(
					self._intrablockbreak.join(subblock)
					for subblock in self._block
				) + "\n").encode()
			)

			# Flush the block/#stored.
			self._block = []
			self._stored = 0

	
	def __enter__(self):
		"""
		Required context management magic method.
		"""
		return self


	def __exit__(self, exc_type, exc_value, exc_tb):
		"""
		Required context management magic method: writes what's left of the cache
		to file, then closes the writer.
		"""
		# If we've finished storing but there's still stuff left in the cache,
		# write that to file too.
		if len(self._block):
			self._writer.write(
				(self._interblockbreak.join(
					self._intrablockbreak.join(subblock)
					for subblock in self._block
				) + "\n").encode()
			)
		
		self._writer.close()


class Player():
	"""
	Safely "replays" recorded `Chain` output.
	"""
	def __init__(self): pass

	def playback(self, fp:str, steps=1000, dtypes=(int, int, int)):
		"""
		Initialize playback; context management.

		Args:
			fp (str): File in which records are stored.
			steps (int=1000): How many steps in the chain? Only relevant for
				displaying the progress bar.
			dtypes (tuple=(int, int, int)): NumPy dtypes for reporting data.

		Returns:
			This Player object.
		"""
		# Configure filepath, reader.
		self._fp = fp
		self._reader = lz4.frame.LZ4FrameFile(f"{self._fp}", mode="rb")
		self._dtypes = dtypes

		# Similar setup to the Recorder.
		self._current = None
		self._currentized = False

		self._loaded = []
		self._remaining = 0
		self._blocksize = 0
		self._configurationsize = -1

		# Compression markers.
		self._intrablockbreak = '<<#>>'
		self._interblockbreak = '<<<##>>>'
		self._interlistbreak = '<>'
		self._listsep = ' '

		# Number of steps (progress bar only).
		self._steps = steps

		# Enter context management.
		return self.__enter__()
	

	def __iter__(self):
		"""
		This `Player` is a generator which yields states of the recorded chain.
		"""
		return self


	def __next__(self):
		"""
		Iteration magic method.
		"""
		# If there are no more configurations remaining in this block, load the
		# next block.
		if self._remaining == 0:
			# Check whether we're EOF.
			try:
				unconfigured = self._reader.readline()
				assert unconfigured != b''
			except:
				raise StopIteration
			
			# Split on block markers...
			unpackedlines = unconfigured.decode().strip().split(self._interblockbreak)

			# ... then on intra-block markers (indicating iterations) ...
			unpackediterations = [l.split(self._intrablockbreak) for l in unpackedlines]

			# ... then on intra-list markers (indicating different arrays within each iteration)...
			self._loaded = [[s.split(self._interlistbreak) for s in l] for l in unpackediterations]
			self._loaded = [
				[
					[np.fromstring(i, sep=self._listsep, dtype=int), np.fromstring(v, sep=self._listsep, dtype=int)]
					for i, v in iteration
				]
				for iteration in self._loaded	
			]

			# Compute the number of configurations that need to be reported and
			# the size of the currently-loaded block (in iterations).
			self._remaining = len(self._loaded)
			self._blocksize = len(self._loaded)

			# Check if we need to load an initial configuration; this only
			# happens on the first step.
			if not self._currentized:
				self._current = tuple([self._loaded[0][i][1] for i in range(len(self._loaded[0]))])
				self._configurationsize = len(self._current)
				self._currentized = True

		# Load the next configuration, modify the current one.
		nextup = self._loaded[self._blocksize-self._remaining]

		for i in range(self._configurationsize):
			indices, values = nextup[i]
			values = values.astype(self._dtypes[i])
			self._current[i][indices] = values
		
		# Decrement the number of configurations remaining in the block, and
		# return the current(ly modified) configuration.
		self._remaining -= 1

		return self._current

	def __enter__(self):
		"""
		Required context management magic method.
		"""
		return self


	def __exit__(self, exc_type, exc_value, exc_tb):
		"""
		Required context management magic method; kills the reader and the file.
		"""
		self._reader.close()

	
	def progress(self):
		from tqdm.auto import tqdm
		return tqdm(self, total=self._steps)

