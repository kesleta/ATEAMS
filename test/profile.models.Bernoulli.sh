#!/bin/zsh

autoload colors; colors
start="${1:-6}"
stop="${2:-9}"
dim="${3:-4}"

echo "_____________________"
echo "| PROFILE BERNOULLI | ➭➭➭ results in profiles/Bernoulli"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %s" "     " "SCALE" "DIM"
echo

for ((L=$start; L<$stop; L++)); do
	python profile.models.Bernoulli.py $dim $L &
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