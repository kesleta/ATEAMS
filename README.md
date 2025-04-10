

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14284172.svg)](https://doi.org/10.5281/zenodo.14284172)

ATEAMS is designed for high-performance simulation of the _Potts_ and _random-cluster_ models in cubical complexes of arbitrary dimension and scale. The linear algebra subroutines supporting these programs are tailored to this application — matrix reduction over finite fields — are implemented in C/C++ using cached 16-bit integer arithmetic. The option to perform these operations in parallel, making computations on large systems manageable with high-performance computers with extremely low memory overhead.

## Installation
Clone this repository by `git clone https://github.com/apizzimenti/ateams.git`, navigate into the project's root, and run `make install` in your favorite terminal. This installs the `ateams` package (and its dependencies) globally in development mode, so any changes you make to the source are reflected system-wide.

### Notes:
* **Read the installation note in `setup.py` to ensure PHAT is correctly installed on your system.**
* **Ensure you have a C/C++ compiler that supports the `-fopenmp` compilation and linking flags.**
* **Before performing sparse/parallel computations, run `make profile` to determine whether sparse/parallel computations make sense for your machine.**

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
@software{ateams,
    title={{ATEAMS: Algebraic Topology-Enabled AlgorithMs for Spin systems}},
    author={Duncan, Paul and Pizzimenti, Anthony E. and Schweinhart, Benjamin},
    url={https://github.com/apizzimenti/ateams},
    version={1.0.2},
    doi={10.5281/zenodo.14284172}
}
```
