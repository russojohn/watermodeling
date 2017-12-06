#!/usr/bin/env bash

init="input/initial_conditions.dat"
home="~/Desktop/project2"
itr=1000000

cd $home
rm -r output/*

for temp in `seq 0.0 0.1 2.0`
do
	folder=temp_${temp}

	cd output
	mkdir $folder
	cd ..
	./build/start $init $temp $itr

	mv pressure.dat output/$folder/
	mv energy.dat output/$folder/
	mv output_* output/$folder/

	awk -vT=$temp '{sum+=$2}END{print T,sum/NR}' output/$folder/energy.dat >> output/avg_energy.txt
done
