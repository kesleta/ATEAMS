
SPARSE=("dense" "sparse")
PARALLEL=("serial" "parallel")
autoload colors; colors

for ((sparse=0; sparse<2; sparse++)); do
	for ((parallel=0; parallel<2; parallel++)); do
		printf "%-25s" "${SPARSE[$(($sparse+1))]}, ${PARALLEL[$(($parallel+1))]}"
		python test.arithmetic.matrices.py $sparse $parallel > profiles/local/${SPARSE[$(($sparse+1))]}.${PARALLEL[$(($parallel+1))]}.txt &
		wait $!
		if [ $? != 0 ]; then
			echo "$fg[red]FAIL$reset_color"
		else
			echo "$fg[green]PASS$reset_color"
		fi
	done
done
