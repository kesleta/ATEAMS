#!/bin/zsh

autoload colors; colors
start="${1:-3}"
stop="${2:-4}"
dim="${3:-4}"
fields=(2 3 5 7)


source .vars

echo "_________________________"
echo "| PROFILE SWENDSEN-WANG | ➭➭➭ results in profiles/SwendsenWang"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %s" "     " "SCALE" "DIM" "FIELD"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do	
		python profile.models.SW.py $dim $L $field &
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