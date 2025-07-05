#!/bin/zsh

autoload colors; colors
start="${1:-9}"
stop="${2:-12}"
cores="${3:-2}"
fields=(2 3 5 7)

source .vars

echo "____________________"
echo "| PROFILE NIENHUIS | ➭➭➭ results in profiles/Nienhuis"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %s" "     " "SCALE" "FIELD"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.NH.py $L $field &
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

echo -e $PASSDESC
echo -e $WARNDESC
echo -e $FAILDESC
echo
echo
echo