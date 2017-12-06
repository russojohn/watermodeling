#!/bin/bash



for T in 0.09 0.1 0.15 0.3 0.45 0.6
do
	mkdir T_$T

	cd T_$T

	cp ../initial_conditions1.dat .

	../../watermodeling/water initial_conditions1.dat $T > output.txt &

	cd ..

done
