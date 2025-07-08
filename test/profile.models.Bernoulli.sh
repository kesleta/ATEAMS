#!/bin/zsh

start="${1:-6}"
stop="${2:-9}"
dim="${3:-4}"

echo "_____________________"
echo "| PROFILE BERNOULLI | ➭➭➭ results in profiles/Bernoulli"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

source .vars

printf "%-5s %-${W}s %-${W}s %s" "" $SCALE $DIM $CELLS
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
echo
echo