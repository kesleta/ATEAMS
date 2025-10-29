
SHELL := /bin/zsh
DIR := $(dir $(abspath $(firstword $(MAKEFILE_LIST))))


###############
#### BUILD ####
###############
clean:
	@rm -rf ./build


quick: FORCE headers
	@python setup.py build_ext --inplace


build: clean PHATMethods LinBoxMethods
	@python setup.py build_ext --inplace > build.log 2>&1


headers:
	@sudo cp -r ateams/common.h /usr/local/include/ATEAMS/
	@sudo cp -r ateams/arithmetic/include/PHAT /usr/local/include/phat
	@sudo cp -r ateams/arithmetic/PHATMethods.h /usr/local/include/ATEAMS/
	@sudo cp -r ateams/arithmetic/LinBoxMethods.h /usr/local/include/ATEAMS/


PHAT_LFLAGS = -I/usr/local/include -shared -fPIC -Wl,-rpath,/usr/local/lib
PHAT_CFLAGS = -O2 -std=c++17

PHATMethods: headers
	@sudo clang++ $(PHAT_CFLAGS) $(PHAT_LFLAGS) -o /usr/local/lib/libPHATMethods.so ateams/arithmetic/PHATMethods.cpp



LINBOX_LFLAGS = -I/usr/local/include -L/usr/local/lib `pkg-config --libs linbox` -lspasm -shared -fPIC -Wl,-rpath,/usr/local/lib
LINBOX_CFLAGS = -O2 -std=c++17

LinBoxMethods: headers
	@sudo clang++ $(LINBOX_CFLAGS) $(LINBOX_LFLAGS) -o /usr/local/lib/libLinBoxMethods.so ateams/arithmetic/LinBoxMethods.cpp 



#################
#### PROFILE ####
#################
Glauber: FORCE
	@cd test && ./profile.models.Glauber.sh 9 12 4
	@cd test && ./profile.models.Glauber.sh 99 102 2

SwendsenWang: FORCE
	@cd test && ./profile.models.SW.sh 4 7 4
	@cd test && ./profile.models.SW.sh 49 52 2

Nienhuis: FORCE
	@cd test && ./profile.models.NH.sh 49 52 2
	@cd test && ./profile.models.NH.sh 9 12 3

InvadedCluster: FORCE
	@cd test && ./profile.models.IC.sh 4 5 4
	@cd test && ./profile.models.IC.sh 19 22 2

Bernoulli: FORCE
	@cd test && ./profile.models.Bernoulli.sh 4 7 4
	@cd test && ./profile.models.Bernoulli.sh 19 22 2


profile: Glauber SwendsenWang Nienhuis InvadedCluster

test: FORCE
	

gauntlet: FORCE test profile

# Kills all the profiles (in case something goes wrong, or is taking too long).
killall: FORCE
	@ps aux | grep -e python -e make -e profile.models | awk '{print $$2}' | xargs kill




#######################
#### DOCUMENTATION ####
#######################
tables: FORCE
	@echo "This will hang if you aren't logged in as anthony on Pangolin."
	@cd test/stats && ./stats.sync.sh && ./stats.tables.sh

docs: FORCE quick
	@pdoc ateams --force --html --template-dir docs/templates --output-dir=docs
	@rsync -a docs/ateams/ docs/
	@rm -rf docs/ateams
	@open "file://$(DIR)docs/index.html"

refs: FORCE
	@pandoc -t markdown_strict --citeproc _refs.md -o refs.md

contribute: build profile docs refs





#################
#### INSTALL ####
#################
install: dependencies _install profile

dependencies: FORCE headers
	@pip install -r requirements.txt

_install: FORCE build
	python setup.py develop

FORCE:
