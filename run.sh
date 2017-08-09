#!/bin/bash

g++-7 -O3 -std=c++14 -fopenmp -Itrng-4.19 -Ltrng-4.19/src/.libs -ltrng4 1D_new_parallel.cpp

type_of_rate_func=0
n=0

#Remove the file if testing.csv file exists
#Otherwise, create the file
#make sure the file is freshly created for each time of execution

if [ -f testing.csv ]; 
	then
		rm testing.csv
		if [ ! -f testing.csv ]; 
		then
			echo "Revomed Successfully!"
		else
			echo "Oops! Something Wrong!"
		fi
fi

while [ $n -le 10 ]
do
	echo "Welcome $n loop"
	./a.out >> testing.csv
	#./a.out
	n=$(( n + 1 ))
done

unset n
unset type_of_rate_func