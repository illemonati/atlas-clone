#! /bin/bash

. $(dirname $0)/find_atlas

mkdir majorMinor
file=majorMinor/times
if [ ! -f "$file" ]; then
	echo -e "majorMinor" > majorMinor/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task majorMinor --glf GLF/GLF.glf.gz --fixedSeed 0 --logFile majorMinor/majorMinor.out --minMAF 0.05 --out majorMinor/majorMinor

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> majorMinor/times
