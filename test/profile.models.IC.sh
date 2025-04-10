#!/bin/zsh

start="${1:-3}"
stop="${2:-4}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-64}"
cores="${5:-2}"

echo "___________________________"
echo "| PROFILE INVADED-CLUSTER | ➭➭➭ results in profiles/InvadedCluster"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-10s %-10s %-5s %-5s %-5s %s" "SIZE" "SPARSE" "PARALLEL" "MIN" "MAX" "CORES"
echo

for ((L=$start; L<$stop; L++)); do
	for ((sparse=0; sparse<2; sparse++)); do
		for ((parallel=0; parallel<2; parallel++)); do
			python profile.models.IC.py $L $sparse $parallel $minBlockSize $maxBlockSize $cores 0 &
			wait $!
		done
	done
done

echo
echo