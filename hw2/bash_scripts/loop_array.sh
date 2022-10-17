#!/bin/bash

# organize FOR loop printing the elements of array
echo "print array:"
arr=("128" "256" "312" "356" "4l2" "5l2" "6l2" "768" "8l2" "9l2")
for i in `seq 1 ${#arr[@]}`; do
    echo ${arr[$i-1]}
done
