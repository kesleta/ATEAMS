

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14284172.svg)](https://doi.org/10.5281/zenodo.14284172)

ATEAMS is designed for high-performance simulation of the _Potts_ and _random-cluster_ models in cubical complexes of arbitrary dimension and scale. The linear algebra subroutines supporting these programs are tailored to this application — matrix reduction over finite fields — and are implemented in C/C++ using cached 16-bit integer arithmetic.

## Installation
Clone this repository by `git clone https://github.com/apizzimenti/ateams.git`, navigate into the project's root, and run `make install` in your favorite terminal. This installs the `ateams` package (and its dependencies) globally in development mode, so any changes you make to the source are reflected system-wide.

### Notes:
* **Read the installation note in `setup.py` to ensure PHAT is correctly installed on your system.**
* **Use the Clang/LLVM compiler family, as some features break when using GCC compilers.**
* **Before performing sparse/parallel computations, run `make profile` to determine whether sparse/parallel computations make sense for your machine.**
* For C++ support, you must install GMP, OpenBLAS (or another BLAS library), Givaro, fflas-ffpack, and LinBox. In general, the auto-installation tools provided by LinBox do not work, so these installations will have to be done manually. It is horrible.

## Example Use

Simulating the $1$-dimensional plaquette random cluster model on a $10 \times 10 \times 10$ cubical $3$-torus with coefficients in the finite field $\mathbb F_3$ looks like:

```python
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain

complex = Cubical().fromCorners([10]*3, field=3)
IC = InvadedCluster(complex, dimension=1)

for (spins, occupied, satisfied) in Chain(IC, steps=10):
    <do whatever>
```

and the $2$-dimensional plaquette Swendsen-Wang algorithm at criticality on a scale $10$ cubical $4$-torus with coefficients in $\mathbb F_5$ looks like

```python
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import critical
from ateams import Chain

complex = Cubical().fromCorners([10]*4, field=5)
SW = SwendsenWang(complex, dimension=2, temperature=critical(complex.field))

for (spins, occupied, satisfied) in Chain(SW, steps=1000):
    <do whatever>
```
The plaquette Swendsen-Wang is, besides Glauber, the most efficient implementation in this library; the above chain (excluding the time required to construct the lattice) runs in ~10 seconds using LinBox on an Apple M2. There are ${\approx}2.4 \times 10^9$ total entries in the boundary matrix for this particular $4$-torus, but only ${\approx}2.4 \times 10^5$ are nonzero, for a density of ${\approx}0.01\%$; the LinBox features _immensely_ reduce the time required to perform the matrix-reduction computations.


You can turn on a progress bar for your simulation using the

```python
...
for (spins, occupied, satisfied) in Chain(HP, steps=10).progress():
    <do whatever>
```
pattern.

## Contributing

* **Do not push directly to this repository: use the pull request model.**
* Follow the standard practices already used in this library, including documentation using [PEP8](https://peps.python.org/pep-0008/) and [PEP257](https://peps.python.org/pep-0257/) guidelines.
* When creating mathematical computation routines, create a testing file in the `test` directory following the `test.<submodule>.<routine>.py`/`test.<submodule>.<routine>.sh` convention. To run existing tests, run `make test`; to run tests you design, add them to the `test` recipe in the `Makefile`. **Please ensure that your routines are tested against ground truth; for example, test new matrix reduction routines against `NumPy`/`SciPy`/`Galois` routines, not against routines already in this library.** _(For examples, take a look in the `test` directory.)_
* When creating new simulation models or new computation routines, create a profiling file in the `test` directory following the `<model-or-routine>.py`/`profile.<submodule>.<model-or-routine>.py`/`profile.<submodule>.<model-or-routine>.sh` convention. To profile existing code, run `make profile`; to run profiles you design, add them to the `profile` recipe in the `Makefile`. _(For examples, take a look in the `test` directory.)_
* To run all tests and all profiles, run `make gauntlet`.
* Before opening a new pull request, run `make contribute` to perform a clean rebuild of the C/C++ backend and documentation.

## Citing

### BibTeX
```bibtex
@software{ATEAMS,
    title={{ATEAMS: Algebraic Topology-Enabled AlgorithMs for Spin systems}},
    author={Duncan, Paul and Pizzimenti, Anthony E. and Schweinhart, Benjamin},
    url={https://github.com/apizzimenti/ATEAMS},
    version={2.0.0},
    doi={10.5281/zenodo.14284172}
}
```
