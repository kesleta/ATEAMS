
contribute: build test profile docs


gauntlet: quick test profile


quick:
	@python setup.py build_ext --inplace


build: clean
	python setup.py build_ext --inplace > build.log 2>&1 


test: FORCE
	@cd test && ./test.arithmetic.matrices.sh
	@cd test && ./test.arithmetic.persistence.sh


profile: FORCE
	@cd test && ./profile.models.NH.sh 3 4 32 64 2
	@cd test && ./profile.models.SW.sh 3 4 32 64 2
	@cd test && ./profile.models.IC.sh 3 4 32 64 2


docs: FORCE quick
	./docs.sh

install: _install gauntlet docs


_install: FORCE build
	pip install -e . --config-settings editable_mode=compat

clean:
	@rm -f ateams/*.c*
	@rm -f ateams/*.o
	@rm -f ateams/*.so
	@rm -f ateams/arithmetic/*.c*
	@rm -f ateams/arithmetic/*.html
	@rm -f ateams/arithmetic/*.o
	@rm -f ateams/arithmetic/*.so
	@rm -rf ./build


FORCE:
