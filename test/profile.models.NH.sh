
start="${1:-3}"
stop="${2:-4}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-64}"
cores="${5:-2}"

echo "____________________"
echo "| PROFILE NIENHUIS |"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-10s %-10s %-10s %-10s %-10s %-10s %s" "SIZE" "SPARSE" "PARALLEL" "MINBLOCK" "MAXBLOCK" "CORES"
echo

for ((L=$start; L<$stop; L++)); do
	for ((sparse=0; sparse<2; sparse++)); do
		for ((parallel=0; parallel<2; parallel++)); do
			python profile.models.NH.py $L $sparse $parallel $minBlockSize $maxBlockSize $cores 0 &
			wait $!
		done
	done
done
