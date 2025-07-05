#!/bin/zsh

autoload colors; colors
start="${1:-6}"
stop="${2:-9}"
dim="${3:-4}"

source .vars

echo "_____________________"
echo "| PROFILE BERNOULLI | ➭➭➭ results in profiles/Bernoulli"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %s" "     " "SCALE" "DIM"
echo

for ((L=$start; L<$stop; L++)); do
	python profile.models.Bernoulli.py $dim $L &
	wait $!
	
	case $? in
		1)
			echo -e "$UP$FAIL";;
		0)
			echo -e "$UP$PASS";;
		*)
			echo -e "$UP$WARN";;
	esac
done
echo

echo -e $PASSDESC
echo -e $WARNDESC
echo -e $FAILDESC
echo
echo
echo