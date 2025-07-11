

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14284172.svg)](https://doi.org/10.5281/zenodo.14284172)
[![docs](https://img.shields.io/badge/%E2%93%98-Documentation-%230099cd)](https://apizzimenti.github.io/ATEAMS/) 

**A**lgebraic **T**opology-**E**nabled **A**lgorith**M**s for **S**pin systems (**ATEAMS**) is a software suite designed for high-performance simulation of generalized _Potts_ and _random-cluster_ models in combinatorial complexes of arbitrary dimension and scale. The linear algebra subroutines supporting these programs are tailored to this application — matrix reduction over finite fields — using [LinBox](https://github.com/linbox-team/linbox) and PHAT [[1](https://www.sciencedirect.com/science/article/pii/S0747717116300098),[2](https://bitbucket.org/phat-code/phat/src/master/)].

[**Install dependencies**](#dependencies) $\longrightarrow$ [**Install ATEAMS**](#installing-ateams) $\longrightarrow$ [**Documentation**](https://apizzimenti.github.io/ATEAMS/) $\longrightarrow$ [**Contributing**](#contributing) $\longrightarrow$ [**Citing**](#citing)

## Demo

Simulating the $1$-dimensional plaquette random cluster model on a $10 \times 10 \times 10$ cubical $3$-torus with coefficients in the finite field $\mathbb F_3$ looks like:

```python
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain

field = 3
C = Cubical().fromCorners([10]*3)
IC = InvadedCluster(C, dimension=1, field=field)

for (spins, occupied, satisfied) in Chain(IC, steps=10):
	pass
```

and the $2$-dimensional plaquette Swendsen-Wang algorithm at criticality on a scale $10$ cubical $4$-torus with coefficients in $\mathbb F_5$ looks like

```python
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import critical
from ateams import Chain

field = 5
C = Cubical().fromCorners([10]*4)
SW = SwendsenWang(C, dimension=2, field=field, temperature=critical(field))

for (spins, occupied) in Chain(SW, steps=1000):
	pass
```
[`SwendsenWang`](https://apizzimenti.github.io/ATEAMS/models/index.html#ateams.models.SwendsenWang) is, after [`Glauber`](https://apizzimenti.github.io/ATEAMS/models/index.html#ateams.models.Glauber), the most efficient implementation in ATEAMS. The above chain terminates in ~19 seconds on an Apple M2. Additional performance information for each model is included in [the documentation](https://apizzimenti.github.io/ATEAMS/models/index.html).


You can turn on a progress bar for your simulation using the

```python
for (spins, occupied, satisfied) in Chain(HP, steps=10).progress():
	pass
```
pattern.

## Installation

### Prerequisites
1. Patience. 
2. Python $\geq$ 3.9. To manage different versions of Python on your machine, we recommend [pyenv](https://github.com/pyenv/pyenv).
2. A C/C++ compiler. **Clang is recommended; please ensure your machine's default compiler is Clang.**
3. [GNU Make](https://www.gnu.org/software/make/) (or a Windows variant...) to build and install the library, and to build any changes you might make. _(This is optional for Windows users. The recipes in the Makefile can be performed manually, but it will take much longer. If that's undesirable, the Linux Subsystem for Windows is a useful workaround.)_
4. Standard tools (pkg-config, autoconf, libtool, etc.) for building and maintaining C++ libraries. For Windows, the Visual Studio BuildTools (which include Clang/LLVM) may be required.
5. If you want to keep your sanity, a computer running macOS or a flavor of Linux.

### Installing ATEAMS

0. [Install all dependencies](#dependencies).
1. Clone this repository.
2. Navigate into the ATEAMS directory.
3. Run `make install`.

In summary,

```
$ git clone https://github.com/apizzimenti/ATEAMS.git
  ...

$ cd ATEAMS
$ make install
```

**Should you run into errors,** the `make install` recipe performs the following operations in the order they're listed:

1. Attempts to compile the Python $\leftrightarrow$ Cython $\leftrightarrow$ LinBox C++ interface at `ATEAMS/ateams/arithmetic/LinBoxMethods.cpp`, building it as a shared library and storing it at `/usr/local/lib/libLinBoxMethods.so`.
2. Attempts to compile the Python $\leftrightarrow$ Cython $\leftrightarrow$ LinBox C++ interface at `ATEAMS/ateams/arithmetic/PHATMethods.cpp`, building it as a shared library and storing it at `/usr/local/lib/libPHATMethods.so`.
3. Attempts to compile the Cython components of ATEAMS, spitting the log into a file called `build.log`.
4. Runs `setup.py` and installs the ATEAMS package as a local development package, so it is importable system-wide.
5. Tests arithmetic and profiles the five main models of the library in varying configurations.

**Done manually, these steps are:**

```
$ rm -rf ./build
$ sudo clang++ `pkg-config --libs linbox` -shared -fPIC -std=c++17 -o /usr/local/lib/libLinBoxMethods.so ateams/arithmetic/LinBoxMethods.cpp -v -O3 -ffast-math
  ...

$ sudo cp -r ateams/arithmetic/include/PHAT /usr/local/include/phat
$ sudo clang++ -shared -fPIC -std=c++17 -o /usr/local/lib/libPHATMethods.so ateams/arithmetic/PHATMethods.cpp -v -O3 -ffast-math
  ...

$ python setup.py build_ext --inplace > build.log 2>&1 
$ python setup.py develop
  ...

$ cd test
$ ./profile.models.Glauber.sh 19 22 4
$ ./profile.models.Glauber.sh 999 1002 2
$ ./profile.models.SW.sh 4 7 4
$ ./profile.models.SW.sh 499 502 2
$ ./profile.models.NH.sh 49 52
$ ./profile.models.IC.sh 4 7 4
$ ./profile.models.IC.sh 19 22 2
$ ./profile.models.Bernoulli.sh 4 7 4
$ ./profile.models.Bernoulli.sh 19 22 2
  ...
```

### Dependencies
ATEAMS relies on LinBox. [This link goes to our GitHub fork of LinBox](https://github.com/apizzimenti/linbox), which addresses a small preconditioning bug and modifies its numerical instability warning system so Cython knows when problems arise; otherwise, the library is unchanged from [its original source](https://github.com/linbox-team/linbox). LinBox relies on [fflas-ffpack](https://github.com/linbox-team/fflas-ffpack), [Givaro](https://github.com/linbox-team/givaro), [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS), and [GMP](https://gmplib.org/). **To get the most out of this toolkit, we highly recommended that you install these dependencies.**

#### GMP
1. [Download GMP 6.3.0 from here](https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz).
2. Follow the [installation instructions here](https://gmplib.org/manual/Installing-GMP), passing the `--enable-cxx` flag to the `./configure` script and setting the install prefix to `/usr/local` (or wherever you'd like). In summary, the installation looks like

```
$ wget -c gmplib.org/download/gmp/gmp-6.3.0.tar.xz -O - | tar -xJ
$ cd gmp-6.3.0
$ ./configure --prefix=/usr/local --enable-cxx
  ...

$ make; make check; sudo make install
$ pkg-config --libs gmp gmpxx
  -L/usr/local/lib -lgmpxx -lgmp
```


#### OpenBLAS
In general, follow the [installation instructions here](http://www.openmathlib.org/OpenBLAS/docs/install/). In particular,

* for macOS users, [the openblas formula on homebrew](https://formulae.brew.sh/formula/openblas) is recommended. It will be installed wherever formulae are typically installed on your computer (e.g. `/opt/homebrew/Cellar/openblas/0.3.29/lib`). _If you choose this option, you are done installing OpenBLAS._
* for Linux (e.g. Ubuntu), it's recommended to build OpenBLAS from source. The latest version known to work with all following dependencies is 0.3.29, [which can be found here](https://github.com/OpenMathLib/OpenBLAS/releases/tag/v0.3.29). After downloading and decompressing the tarball, navigate into the directory and run `make`.

Regardless of how you install, OpenBLAS should register with `pkg-config`. In summary: on macOS,
```
$ brew install openblas
$ pkg-config --libs openblas
  -L/opt/homebrew/Cellar/openblas/0.3.29/lib -lopenblas
```

or on Linux,

```
$ wget -c github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29.tar.gz -O - | tar -xz
$ cd OpenBLAS-0.3.29
$ make
  ...

$ pkg-config --libs openblas
  -L/usr/local/lib -lopenblas
```


#### Givaro

To stave off bugs, we recommend cloning [the current main branch of Givaro](https://github.com/linbox-team/givaro) and building from source.

1. Clone the Givaro repository.
2. Navigate into the Givaro directory and run the `autogen` script, passing your preferred install prefix as argument.
3. Run `make` and `sudo make install`.

In summary,

```
$ git clone https://github.com/linbox-team/givaro.git
  Cloning into 'givaro'... done.

$ cd givaro
$ ./autogen.sh --prefix=/usr/local
  ...

$ make; sudo make install
$ pkg-config --libs givaro
  -L/usr/local/lib -lgivaro -lgmpxx -lgmp
```

#### fflas-ffpack

This package can be slightly trickier, and may need some convincing that the previous dependencies actually exist on your system. As with Givaro, we recommend cloning [the current main branch of fflas-ffpack](https://github.com/linbox-team/fflas-ffpack) and building from source.

1. Clone the fflas-ffpack repository.
2. Navigate into the fflas-ffpack directory and run the `autogen` script with the following options:
	* Ubuntu, Debian, Mint, etc.
		* `--prefix=/usr/local`
		* ``--with-blas-libs="`pkg-config --libs openblas`"``
		* ``--with-blas-cflags="`pkg-config --cflags openblas`"``
	* macOS:
		* `--prefix=/usr/local`
		* ``--with-blas-libs="-framework Accelerate"``
		* (it's possible to use the `pkg-config` arguments to ``--with-blas-libs`` and ``--with-blas-cflags`` here, but the LinBox team recommends using the Accelerate framework.)
3. Run `make`.
4. _Optionally_, run `make autotune`.
5. Run `sudo make install; make check`.
6. If/when things go wrong, check [the fflas-ffpack installation notes](https://github.com/linbox-team/fflas-ffpack/blob/master/INSTALL).

In summary,

```
$ git clone https://github.com/linbox-team/fflas-ffpack.git
  Cloning into 'fflas-ffpack'... done.

$ cd fflas-ffpack
$ ./autogen.sh --prefix=/usr/local --with-blas-libs=<libs> --with-blas-cflags=<cflags>
  ...

$ make; make autotune
  ...

$ sudo make install; make check
  ...

$ pkg-config --libs fflas-ffpack
  -L/<path-to-openblas> -lopenblas -L/usr/local/lib -lgivaro -lgmpxx -lgmp
```


#### LinBox

The end is in sight! Here, we recommend cloning [the `bug/bad-checks` branch of our forked LinBox repository](https://github.com/apizzimenti/linbox.git) and building from source.

1. Clone LinBox.
2. Navigate into the LinBox directory and run the `autogen` script specifying the install prefix.
3. Run `make` and `sudo make install`.

In summary,

```
$ git clone https://github.com/apizzimenti/linbox.git
  Cloning into 'linbox'... done.

$ cd linbox
$ git checkout bug/bad-checks
$ ./autogen.sh --prefix=/usr/local
  ...

$ make; sudo make install
  ...

$ pkg-config --libs linbox
  -L/usr/local/lib -llinbox -L/<path-to-openblas> -lopenblas -lgivaro -lgmpxx -lgmp
```


#### PHAT

We use the [Persistent Homology Algorithms Toolbox (PHAT)](https://bitbucket.org/phat-code/phat) to compute persistence over $\mathbb Z/2\mathbb Z$. Included in the `ATEAMS/ateams/arithmetic/include/PHAT` folder are all the header files for the PHAT library (as of writing) which will be copied to `/usr/local/include/phat` whenever `make install` (or `make PHATMethods`) is executed. We build an additional C++ interface in `PHATMethods.cpp`, which is linked against by the Cython compiler and made available to the Python modules in the library --- specifically, native PHAT lets us remarkably speed up persistence computation in the `Bernoulli` and `InvadedCluster` models.

Now, you can move on to [installing ATEAMS itself](#installing-ateams)!

## Contributing

* **Do not push directly to this repository: use the fork + pull request model.**
* Follow the standard practices already used in this library, including documentation according to [PEP8](https://peps.python.org/pep-0008/) and [PEP257](https://peps.python.org/pep-0257/) guidelines.
* When writing new Cython/C/C++ code, please include its compilation in the `build` recipe of the Makefile. If you want C/C++ code to interface with Cython, ensure `setup.py` is correctly configured to find shared libraries.
* When creating mathematical computation routines, create a testing file in the `test` directory following the `test.<submodule>.<routine>.py`/`test.<submodule>.<routine>.sh` convention. To run existing tests, run `make test`; to run tests you design, add them to the `test` recipe in the `Makefile`. **Please ensure that your routines are tested against ground truth; for example, test new matrix reduction routines against `NumPy`/`SciPy` routines, not against routines already in this library.** _(For examples, take a look in the `test` directory.)_
* When creating new simulation models or new computation routines, create a profiling file in the `test` directory following the `<model-or-routine>.py`/`profile.<submodule>.<model-or-routine>.py`/`profile.<submodule>.<model-or-routine>.sh` convention. To profile existing code, run `make profile`; to run profiles you design, add them to the `profile` recipe in the `Makefile`. _(For examples, take a look in the `test` directory.)_
* To run all tests and all profiles, run `make gauntlet`.
* Before opening a new pull request, run `make contribute` to perform a clean rebuild of the Cython/C/C++ backend and documentation.

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

## References
1. Bauer, Ulrich, Michael Kerber, Jan Reininghaus, and Hubert Wagner. 2017. “PHAT — Persistent Homology Algorithms Toolbox.” *Journal of Symbolic Computation* 78 (January): 76–90. <https://doi.org/10.1016/j.jsc.2016.03.008>, <https://bitbucket.org/phat-code/>.

4. Chen, Chao, and Michael Kerber. 2011. “Persistent Homology Computation with a Twist.” In *Proceedings 27th European Workshop on Computational Geometry*, 11:197–200.

7. Duncan, Paul, and Benjamin Schweinhart. 2025. “Topological Phases in the Plaquette Random-Cluster Model and Potts Lattice Gauge Theory.” *Communications in Mathematical Physics*. [arxiv.org/abs/10.48550/arXiv.2207.08339](https://arxiv.org/abs/10.48550/arXiv.2207.08339).

8. Edelsbrunner, Herbert, and John Harer. 2010. *Computational Topology: An Introduction*. American Mathematical Soc.

8. The LinBox group. 2021. *LinBox*. v1.7.0 ed. <http://github.com/linbox-team/linbox>.

9. Hatcher, Allen. 2002. *Algebraic Topology*. Cambridge University Press.

10. Machta, J., Y. S. Choi, A. Lucke, T. Schweizer, and L. M. Chayes. 1996. “Invaded Cluster Algorithm for Potts Models.” *Physical Review E* 54 (2). <https://doi.org/10.1103/PhysRevE.54.1332>.

11. Swendsen, Robert H., and Jian-Sheng Wang. 1987. “Nonuniversal Critical Dynamics in Monte Carlo Simulations.” *Physical Review Letters* 58 (2). <https://doi.org/10.1103/PhysRevLett.58.86>.

12. Zomorodian, Afra J. 2005. *Topology for Computing*. Cambridge University Press.

