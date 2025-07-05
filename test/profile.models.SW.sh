#!/bin/zsh

autoload colors; colors
start="${1:-3}"
stop="${2:-4}"
dim="${3:-4}"
fields=(2 3 5 7)

UP="\033[F"
DOWN="\033[B"
CLEAR="\033[2K"

echo "_________________________"
echo "| PROFILE SWENDSEN-WANG | ➭➭➭ results in profiles/SwendsenWang"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %s" "     " "SCALE" "DIM" "FIELD"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.SW.py $dim $L $field &
		wait $!
		
		if [ $? != 0 ]; then
			echo -e "$UP$UP$UP$UP$UP$fg[red]FAIL$reset_color"
			echo -e "$DOWN$CLEAR$DOWN$CLEAR$DOWN$CLEAR$UP$UP$UP$UP"
		else
			echo -e "$UP$fg[green]PASS$reset_color"
		fi
	done
	echo
done

echo -e "A test $fg[red]FAIL$reset_color typically results from small system size, small field characteristic, or a poorly conditioned matrix."

echo
echo
echo