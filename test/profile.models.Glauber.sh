#!/bin/zsh

autoload colors; colors
start="${1:-99}"
stop="${2:-102}"
dim="${3:-2}"
fields=(2 3 5 7)

echo "___________________"
echo "| PROFILE GLAUBER | ➭➭➭ results in profiles/Glauber"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %-8s %s" "     " "SCALE" "DIM" "FIELD"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.Glauber.py $field $dim $L &
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