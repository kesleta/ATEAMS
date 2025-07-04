#!/bin/zsh

autoload colors; colors
start="${1:-9}"
stop="${2:-12}"
cores="${3:-2}"

echo "____________________"
echo "| PROFILE NIENHUIS | ➭➭➭ results in profiles/Nienhuis"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %-8s %-8s %s" "     " "SCALE" "PARALLEL" "CORES" "LINBOX"
echo

for ((L=$start; L<$stop; L++)); do
	for ((parallel=0; parallel<2; parallel++)); do
		for ((LinBox=0; LinBox<2; LinBox++)); do
			python profile.models.NH.py $LinBox $L $(($parallel & $(($LinBox < 1)))) $cores 0 &
			wait $!
			
			if [ $? != 0 ]; then
				echo -e "\033[F$fg[red]FAIL$reset_color"
			else
				echo -e "\033[F$fg[green]PASS$reset_color"
			fi
		done
	done
	echo
done

echo
echo
echo