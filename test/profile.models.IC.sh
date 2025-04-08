#!/bin/bash

start="${1:-3}"
stop="${2:-4}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-64}"
cores="${5:-2}"

echo "___________________________"
echo "| PROFILE INVADED-CLUSTER |"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-10s %-10s %-10s %-10s %-10s %-10s %s" "SIZE" "SPARSE" "PARALLEL" "MINBLOCK" "MAXBLOCK" "CORES"
echo

for ((L=$start; L<$stop; L++)); do
	for ((sparse=0; sparse<2; sparse++)); do
		for ((parallel=0; parallel<2; parallel++)); do
			python profile.models.IC.py $L $sparse $parallel $minBlockSize $maxBlockSize $cores 0 &
			wait $!
		done
	done
done

# echo "__________"
# echo "| KERNEL |"
# echo "‾‾‾‾‾‾‾‾‾‾"

# for ((sparse=0; sparse<2; sparse++)); do
# 	for ((parallel=0; parallel<2; parallel++)); do
# 		python test.arithmetic.matrices.py $sparse $parallel > profiles/kernel/${SPARSE[$(($sparse+1))]}.${PARALLEL[$(($parallel+1))]}.txt &
# 		wait $!
# 		if [ $? != 0 ]; then
# 			echo -e "\033[F$fg[red]FAIL$reset_color"
# 		else
# 			echo -e "\033[F$fg[green]PASS$reset_color"
# 		fi
# 	done
# done
