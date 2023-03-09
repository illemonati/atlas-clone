#! /bin/bash

. $(dirname $0)/find_atlas

mkdir recal
file=recal/times
if [ ! -f "$file" ]; then
	echo -e "recal" > recal/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task recal --bam simulate/bamFile.bam --recal "intercept[0.1];quality:polynomial[0.8,0.2,0.01]" --onlyLL --fixedSeed 0 --logFile recal/onlyLL.out --out recal/recal

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> recal/times
