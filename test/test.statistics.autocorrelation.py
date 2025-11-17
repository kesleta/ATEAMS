
import numpy as np
from time import time
from ateams.statistics.autocorrelation import unnormalized, normalized, integrated

X = np.random.uniform(0, 1, size=100000)
print(integrated(X))

