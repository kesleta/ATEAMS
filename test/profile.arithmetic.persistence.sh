#!/bin/bash

start="${1:-3}"
stop="${2:-4}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-64}"
cores="${5:-2}"

for ((L=$start; L<$stop; L++)); do
	for ((sparse=0; sparse<2; sparse++)); do
		for ((parallel=0; parallel<2; parallel++)); do
			echo
			echo
			echo -e "executing complex $L"
			echo -e "\tsparse? \t$sparse"
			echo -e "\tparallel? \t$parallel"
			echo -e "\tminBlockSize \t$minBlockSize"
			echo -e "\tmaxBlockSize \t$maxBlockSize"
			echo -e "\tcores \t\t$cores"
			python profile.arithmetic.persistence.py $L $sparse $parallel $minBlockSize $maxBlockSize $cores &
			wait $!
		done
	done
done
