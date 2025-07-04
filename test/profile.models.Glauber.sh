#!/bin/zsh

autoload colors; colors
start="${1:-9}"
stop="${2:-12}"

echo "___________________"
echo "| PROFILE GLAUBER | ➭➭➭ results in profiles/Glauber"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %s" "     " "SCALE"
echo

for ((L=$start; L<$stop; L++)); do
	python profile.models.Glauber.py $L &
	wait $!
	
	if [ $? != 0 ]; then
		echo -e "\033[F$fg[red]FAIL$reset_color"
	else
		echo -e "\033[F$fg[green]PASS$reset_color"
	fi
done

echo
echo
echo