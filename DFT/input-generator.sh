#!/bin/bash

echo "How many random numbers to generate?"
read maxnum
echo Generating $maxnum integers
for (( c=1; c <= $maxnum; c++ ))
do
    ((number = $RANDOM - 16384))
    echo $number
done

