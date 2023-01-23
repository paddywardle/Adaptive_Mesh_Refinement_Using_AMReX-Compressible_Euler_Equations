#!/bin/bash

file="../data/results.dat"

for i in {2..4}
do
	echo "plot '$file' using 1:$i with linespoints" | gnuplot -persist
done	
