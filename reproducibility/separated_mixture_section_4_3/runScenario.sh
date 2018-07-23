#!/bin/bash

for i in `seq 1 2`;
do
    for j in `seq 1 2`;
    do
	for k in `seq 1 2`;
	do
	    for l in `seq 1 2`;
	    do
		/bin/bash /path/to/runSimulation.sh $i $j $k $l
	    done
	done
    done
done 
