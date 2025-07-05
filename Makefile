
contribute: build test profile docs


gauntlet: quick test profile


quick: fast
	@python setup.py build_ext --inplace


build: clean fast
	python setup.py build_ext --inplace > build.log 2>&1 


test: FORCE
	@cd test && ./test.arithmetic.matrices.sh
	@cd test && ./test.arithmetic.bettis.sh


PHAT: FORCE
	@pip install pybind11 wheel setuptools
	@pip install --use-deprecated=legacy-resolver --no-binary :all: phat


fast:
	@sudo clang++ `pkg-config --libs linbox` -shared -fPIC -std=c++17 -o /usr/local/lib/libLinBoxMethods.so ateams/arithmetic/LinBoxMethods.cpp -v -O3 -ffast-math

# Individual model profiles.

Glauber: FORCE
	@cd test && ./profile.models.Glauber.sh 11 14 4

SwendsenWang: FORCE
	@cd test && ./profile.models.SW.sh 9 12 2

Nienhuis: FORCE
	@cd test && ./profile.models.NH.sh 9 12 2

InvadedCluster: FORCE
	@cd test && ./profile.models.IC.sh 3 5 2

Bernoulli: FORCE
	@cd test && ./profile.models.Bernoulli.sh 3 5 2

profile: Glauber SwendsenWang Nienhuis InvadedCluster Bernoulli



# Build docs.
docs: FORCE quick
	./docs.sh


# Installation recipes.
install: _install gauntlet docs

_install: FORCE build PHAT
	python setup.py develop

clean:
	@rm -rf ./build

FORCE:
