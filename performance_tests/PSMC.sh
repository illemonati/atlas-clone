#! /bin/bash

. $(dirname $0)/find_atlas

mkdir PSMC
file=PSMC/times
if [ ! -f "$file" ]; then
	echo -e "PSMC" > PSMC/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task PSMC --bam simulate/bamFile.bam --fixedSeed 0 --logFile PSMC/PSMC.out --out PSMC/PSMC

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> PSMC/times
