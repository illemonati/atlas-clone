#! /bin/bash

. $(dirname $0)/find_atlas

mkdir mergeRG
file=mergeRG/times
if [ ! -f "$file" ]; then
	echo -e "mergeRG" > mergeRG/times
fi

timeFor10Runs=0

echo "SimReadGroup1 SimReadGroup2" > mergeRG/rgs.txt
echo "SimReadGroup3 SimReadGroup4 SimReadGroup5 SimReadGroup6 SimReadGroup7 SimReadGroup8 SimReadGroup9 SimReadGroup10" >> mergeRG/rgs.txt
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task mergeRG --bam simulate/bamFile.bam --readGroups mergeRG/rgs.txt --logFile mergeRG/mergeRG.out --fixedSeed 0 --out mergeRG/mergeRG

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> mergeRG/times
