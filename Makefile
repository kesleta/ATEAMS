
test: FORCE
	@cd test && zsh test.arithmetic.matrices.sh

quick:
	@python setup.py build_ext --inplace

sparse: quick
	@cd test && zsh test.arithmetic.matrices.sh


build: clean
	python setup.py build_ext --inplace

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