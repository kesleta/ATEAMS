
import numpy as np

MINT = np.int64
"""Global "machine integer" data type; this is equivalent to the C `int` data type."""

SFINT = np.int16
"""\"Short finite-field integer\" data type; used to reduce memory overhead, since we need far fewer bits to sore our data."""

FINT = np.int16
"""\"Finite-field integer\" data type; uses 1/4th the space of a machine integer."""
