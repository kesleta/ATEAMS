

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14284172.svg)](https://doi.org/10.5281/zenodo.14284172)

ATEAMS is designed for high-performance simulation of the _Potts_ and _random-cluster_ models in cubical complexes of arbitrary dimension and scale. The linear algebra subroutines supporting these programs are tailored to this application — matrix reduction over finite fields — are implemented in C/C++ using cached 16-bit integer arithmetic. The option to perform these operations in parallel, making computations on large systems manageable with high-performance computers with extremely low memory overhead.

## Installation
Clone this repository by `git clone https://github.com/apizzimenti/ateams.git`, navigate into the project's root, and run `make install` in your favorite terminal. This installs the `ateams` package (and its dependencies) globally in development mode, so any changes you make to the source are reflected system-wide.

### Notes:
* **Read the installation note in `setup.py` to ensure PHAT is correctly installed on your system.**
* **Ensure you have a C/C++ compiler that supports the `-fopenmp` compilation and linking flags.**

## Example Use

Simulating the $1$-dimensional plaquette random cluster model on a $ 10 \times 10 \times 10 $ cubical $3$-torus with coefficients in the finite field $\mathbb F_3$ looks like:

```python
from ateams.structures import Lattice
from ateams.models import InvadedCluster
from ateams import Chain

L = Lattice().fromCorners([10,10,10], field=3)
HP = InvadedCluster(L, homology=1)

for (spins, occupied, satisfied) in Chain(HP, steps=10):
    <do whatever>
```

You can turn on a progress bar for your simulation using the

```python
...
for (spins, occupied, satisfied) in Chain(HP, steps=10).progress():
    <do whatever>
```
pattern.

## Citing

### BibTeX
```bibtex
@software{ateams,
    title={{ATEAMS: Algebraic Topology-Enabled AlgorithMs for Spin systems}},
    author={Duncan, Paul and Pizzimenti, Anthony E. and Schweinhart, Benjamin},
    url={https://github.com/apizzimenti/ateams},
    version={1.0.0},
    doi={10.5281/zenodo.14284172}
}
```
