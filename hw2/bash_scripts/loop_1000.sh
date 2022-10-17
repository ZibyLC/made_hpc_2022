#!/bin/bash

#  organize FOR loop printing the even numbers only from 100 to 1000
start=100
end=1000
echo "all even numbers from $start to $end"
for n in `seq $start 2 $end`; do
  echo -n "$n "
done
echo ""
