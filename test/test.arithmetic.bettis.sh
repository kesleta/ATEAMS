#!/bin/zsh

autoload colors; colors

echo "_______________"
echo "| TEST BETTIS | ➭➭➭ results in profiles/bettis"
echo "‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾"

for ((sparse=0; sparse<1; sparse++)); do
	python test.arithmetic.bettis.py > profiles/bettis/sparse.serial.txt &
	wait $!
	if [ $? != 0 ]; then
		echo -e "\033[F$fg[red]FAIL$reset_color"
	else
		echo -e "\033[F$fg[green]PASS$reset_color"
	fi
done

echo
echo
