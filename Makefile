
quick:
	@python setup.py build_ext --inplace

build: clean
	python setup.py build_ext --inplace > build.log 2>&1 

test: FORCE
	@cd test && ./test.arithmetic.matrices.sh
	@echo
	@echo
	@cd test && ./test.arithmetic.persistence.sh

profile: FORCE
	@cd test && ./profile.models.IC.sh 3 7 32 64 18
	# @cd test && zsh profile.arithmetic.persistence.sh 5 6 6 32 64 2

sparse: quick
	@cd test && zsh test.arithmetic.matrices.sh

docs: FORCE quick
	./docs.sh

install: FORCE build
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