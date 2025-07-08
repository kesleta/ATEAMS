#!/bin/zsh

start="${1:-99}"
stop="${2:-102}"
dim="${3:-2}"
fields=(2 3 5 7)

echo "___________________"
echo "| PROFILE GLAUBER | ➭➭➭ results in profiles/Glauber"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

source .vars

printf "%-5s %-${W}s %-${W}s %-${W}s %s" "" $SCALE $DIM $FIELD $CELLS
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.Glauber.py $field $dim $L &
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
done

echo
echo