#!/bin/zsh

start="${1:-9}"
stop="${2:-12}"
dim="${3:-2}"
fields=(2 3 5 7)

echo "____________________"
echo "| PROFILE NIENHUIS | ➭➭➭ results in profiles/Nienhuis"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

source .vars

printf "%5s %-${W}s %-${W}s %-${W}s %s" "" $SCALE $DIM $FIELD $CELLS
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.NH.py $dim $L $field &
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