#!/bin/zsh

start="${1:-99}"
stop="${2:-102}"
dim="${3:-2}"
fields=(2 3 5 7)

source .vars


echo "___________________"
echo "| PROFILE GLAUBER | ➭➭➭ results in profiles/Glauber"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %-8s %s" "     " "SCALE" "DIM" "FIELD"
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

echo -e $PASSDESC
echo -e $WARNDESC
echo -e $FAILDESC
echo
echo
echo