#!/bin/zsh


SPARSE=("dense" "sparse")
PARALLEL=("serial" "parallel")
autoload colors; colors

echo "_______________"
echo "| TEST KERNEL | ➭➭➭ results in profiles/kernel"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

for ((sparse=0; sparse<2; sparse++)); do
	for ((parallel=0; parallel<2; parallel++)); do
		python test.arithmetic.matrices.py $sparse $parallel > profiles/kernel/${SPARSE[$(($sparse+1))]}.${PARALLEL[$(($parallel+1))]}.txt &
		wait $!
		if [ $? != 0 ]; then
			echo -e "\033[F$fg[red]FAIL$reset_color"
		else
			echo -e "\033[F$fg[green]PASS$reset_color"
		fi
	done
done

echo
echo
