#!/bin/bash

# compute 100 + 0.5 in bash

A=100
B=0.5
echo "scale=10; $(($A + $B))" | bc -l
