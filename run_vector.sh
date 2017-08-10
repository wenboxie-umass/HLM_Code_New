#!/bin/bash

g++-7 -O3 -std=c++14 -fopenmp -Itrng-4.19 -Ltrng-4.19/src/.libs -ltrng4 1D_new_vector_parallel.cpp -o 1D_new_vector_parallel

if [ -f running_log_vector.txt ];
	then
		rm running_log_vector.txt;
		echo "Removed Successfully"
fi 

n=1

while [ $n -le 2 ]
do
	echo "Iteration: $n"
	echo "Iteration: $n" >> running_log_vector.txt
	./1D_new_vector_parallel $n >> running_log_vector.txt
	n=$(( n + 1 ))
done