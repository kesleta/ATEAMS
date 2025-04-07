
quick:
	@python setup.py build_ext --inplace

build: clean
	python setup.py build_ext --inplace

test: FORCE
	# @cd test && zsh test.arithmetic.matrices.sh
	@cd test && zsh test.arithmetic.persistence.sh

profile: FORCE
	@cd test && zsh profile.models.IC.sh 4 5 32 64 2
	# @cd test && zsh profile.arithmetic.persistence.sh 5 6 6 32 64 2

sparse: quick
	@cd test && zsh test.arithmetic.matrices.sh

docs:
	sh docs.sh


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