
SHELL := /bin/zsh

contribute: build test profile docs


gauntlet: quick test profile


quick: PHATMethods LinBoxMethods
	@python setup.py build_ext --inplace


build: clean PHATMethods LinBoxMethods
	python setup.py build_ext --inplace > build.log 2>&1 


test: FORCE
	@cd test && ./test.arithmetic.matrices.sh
	@cd test && ./test.arithmetic.bettis.sh


PHATMethods:
	@sudo cp -r ateams/arithmetic/include/PHAT /usr/local/include/phat
	@sudo clang++ -shared -fPIC -std=c++17 -o /usr/local/lib/libPHATMethods.so ateams/arithmetic/PHATMethods.cpp -v -O3 -ffast-math


LinBoxMethods:
	@sudo clang++ `pkg-config --libs linbox` -shared -fPIC -std=c++17 -o /usr/local/lib/libLinBoxMethods.so ateams/arithmetic/LinBoxMethods.cpp -v -O3 -ffast-math

# Individual model profiles.

LinBox:
	@cd ../libraries/linbox && make; sudo make install

Glauber: FORCE
	# @cd test && ./profile.models.Glauber.sh 19 22 4
	# @cd test && ./profile.models.Glauber.sh 999 1002 2
	@cd test && ./profile.models.Glauber.sh

SwendsenWang: FORCE
	# @cd test && ./profile.models.SW.sh 4 15 4
	# @cd test && ./profile.models.SW.sh 499 502 2
	@cd test && ./profile.models.SW.sh

Nienhuis: FORCE
	# @cd test && ./profile.models.NH.sh 49 52
	@cd test && ./profile.models.NH.sh

InvadedCluster: FORCE
	# @cd test && ./profile.models.IC.sh 4 8 4
	# @cd test && ./profile.models.IC.sh 19 22 2
	@cd test && ./profile.models.IC.sh

Bernoulli: FORCE
	# @cd test && ./profile.models.Bernoulli.sh 4 9 4
	# @cd test && ./profile.models.Bernoulli.sh 49 52 2
	@cd test && ./profile.models.Bernoulli.sh

profile: Glauber SwendsenWang Nienhuis InvadedCluster Bernoulli

killall: FORCE
	@ps aux | grep -e python -e make -e profile.models | awk '{print $$2}' | xargs kill



# Build docs.

tables: FORCE
	@cd test/stats && ./stats.sync.sh && ./stats.tables.sh

docs: FORCE tables
	./docs.sh


# Installation recipes.
install: _install gauntlet docs

_install: FORCE build
	python setup.py develop

clean:
	@rm -rf ./build

FORCE:
