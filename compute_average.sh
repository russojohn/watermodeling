#!/bin/bash



for T in 0.09 0.1 0.15 0.3 0.45 0.6
do

	cd T_$T

	awk -vT=$T '{sum+=$2}END{print T,sum/NR}' energy.dat

	cd ..

done
