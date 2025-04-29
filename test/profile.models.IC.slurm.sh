
start="${1:-3}"
stop="${2:-8}"
minBlockSize="${3:-32}"
maxBlockSize="${4:-128}"
cores="${5:-12}"

printf "%-10s %-10s %-10s %s" "SPARSE" "PARALLEL" "SIZE" "JOB"
echo

for ((L=$start; L<$stop; L++)); do
	for ((sparse=0; sparse<2; sparse++)); do
		for ((parallel=0; parallel<2; parallel++)); do
			last=$(sbatch -J "pr.$sparse.$parallel.$L" profile.models.IC.slurm $L $sparse $parallel $minBlockSize $maxBlockSize $cores 1 | awk '{print $NF}')
			printf "%-10s %-10s %-10s %s" "$sparse" "$parallel" "$L" "$last"
			echo
		done
	done
done
