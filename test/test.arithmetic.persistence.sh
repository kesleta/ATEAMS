#!/bin/zsh

SPARSE=("dense" "sparse")
PARALLEL=("serial" "parallel")
autoload colors; colors

echo "_______________"
echo "| PERSISTENCE |"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

ABOVE="\033[F"

for ((sparse=0; sparse<2; sparse++)); do
	for ((parallel=0; parallel<2; parallel++)); do
		python test.arithmetic.persistence.py $sparse $parallel > profiles/persistence/${SPARSE[$(($sparse+1))]}.${PARALLEL[$(($parallel+1))]}.txt &
		wait $!
		if [ $? != 0 ]; then
			echo -e "$ABOVE$fg[red]FAIL$reset_color"
		else
			echo -e "$ABOVE$fg[green]PASS$reset_color"
		fi
	done
done
