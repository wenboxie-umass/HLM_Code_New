#!/bin/bash

g++-7 -O3 -std=c++14 -fopenmp -Itrng-4.19 -Ltrng-4.19/src/.libs -ltrng4 1D_new_parallel.cpp

n=0

while [ $n -le 20 ]
do
	echo "Welcome $n loop"
	./a.out
	n=$(( n + 1 ))
done