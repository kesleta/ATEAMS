#!/bin/zsh

start="${1:-3}"
stop="${2:-4}"
dim="${3:-4}"
fields=(2 3 5 7)

source .vars

echo "___________________________"
echo "| PROFILE INVADED-CLUSTER | ➭➭➭ results in profiles/InvadedCluster"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

printf "%-5s %-8s %-8s %-8s %s" "     " "SCALE" "DIM" "FIELD" "PHAT"
echo

for ((L=$start; L<$stop; L++)); do
	for field in ${fields[@]}; do
		python profile.models.IC.py $dim $L $field &
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

width=$(tput cols)
indent=2

echo -e $PASSDESC
echo -e $WARNDESC
echo -e $FAILDESC
echo
echo
echo