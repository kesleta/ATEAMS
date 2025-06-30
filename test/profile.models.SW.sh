#!/bin/zsh

start="${1:-3}"
stop="${2:-4}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-64}"
cores="${5:-2}"

echo "_________________________"
echo "| PROFILE SWENDSEN-WANG | ➭➭➭ results in profiles/SwendsenWang"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-10s %-10s %-5s %-5s %-5s %-10s %s" "SCALE" "SPARSE" "PARALLEL" "MIN" "MAX" "CORES" "LINBOX"
echo

# for ((L=$start; L<$stop; L++)); do
# 	for ((sparse=0; sparse<2; sparse++)); do
# 		for ((parallel=0; parallel<2; parallel++)); do
# 			for ((LinBox=0; LinBox<2; LinBox++)); do
# 				python profile.models.SW.py $LinBox $L $sparse $parallel $minBlockSize $maxBlockSize $cores 0 &
# 				wait $!
# 			done
# 		done
# 	done
# done

for ((L=$start; L<$stop; L++)); do
	for ((LinBox=0; LinBox<1; LinBox++)); do
		python profile.models.SW.py 1 $L 1 0 32 64 2 0 &
		wait $!
	done
done

echo
echo
