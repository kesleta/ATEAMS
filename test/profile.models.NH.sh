#!/bin/zsh

autoload colors; colors
start="${1:-9}"
stop="${2:-12}"
cores="${3:-2}"
fields=(2 3 5 7)

echo "____________________"
echo "| PROFILE NIENHUIS | ➭➭➭ results in profiles/Nienhuis"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %s" "     " "SCALE" "FIELD"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.NH.py $L $field &
		wait $!
		
		if [ $? != 0 ]; then
			echo -e "\033[F$fg[red]FAIL$reset_color"
		else
			echo -e "\033[F$fg[green]PASS$reset_color"
		fi
	done
	echo
done

echo
echo
echo