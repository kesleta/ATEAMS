
contribute: build test profile docs


gauntlet: quick test profile


quick: fast
	@python setup.py build_ext --inplace


build: clean fast
	python setup.py build_ext --inplace > build.log 2>&1 


test: FORCE
	@cd test && ./test.arithmetic.matrices.sh
	@cd test && ./test.arithmetic.bettis.sh


PHAT:
	@sudo cp -r ateams/arithmetic/include/PHAT /usr/local/include/phat
	@sudo clang++ -shared -fPIC -std=c++17 -o /usr/local/lib/libPHATMethods.so ateams/arithmetic/PHATMethods.cpp -v -O3 -ffast-math


fast: PHAT
	@sudo clang++ `pkg-config --libs linbox` -shared -fPIC -std=c++17 -o /usr/local/lib/libLinBoxMethods.so ateams/arithmetic/LinBoxMethods.cpp -v -O3 -ffast-math

# Individual model profiles.

Glauber: FORCE
	@cd test && ./profile.models.Glauber.sh 11 12 4
	@cd test && ./profile.models.Glauber.sh

SwendsenWang: FORCE
	@cd test && ./profile.models.SW.sh 12 16 4
	@cd test && ./profile.models.SW.sh 499 502 2

Nienhuis: FORCE
	@cd test && ./profile.models.NH.sh 99 102

InvadedCluster: FORCE
	@cd test && ./profile.models.IC.sh 4 5 4
	@cd test && ./profile.models.IC.sh 69 72 2

Bernoulli: FORCE
	@cd test && ./profile.models.Bernoulli.sh 6 7 4
	@cd test && ./profile.models.Bernoulli.sh 9 11 2

profile: Glauber SwendsenWang Nienhuis InvadedCluster Bernoulli



# Build docs.
docs: FORCE quick
	./docs.sh


# Installation recipes.
install: _install gauntlet docs

_install: FORCE build
	python setup.py develop

clean:
	@rm -rf ./build

FORCE:
