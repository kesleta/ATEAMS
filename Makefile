
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


profile: FORCE
	@cd test && ./profile.models.NH.sh 9 12 2
	@cd test && ./profile.models.SW.sh 3 5 2
	@cd test && ./profile.models.IC.sh 3 5 2
	@cd test && ./profile.models.Glauber.sh 9 12


docs: FORCE quick
	./docs.sh

install: _install gauntlet docs


_install: FORCE build PHAT
	python setup.py develop

clean:
	@rm -rf ./build


FORCE:
