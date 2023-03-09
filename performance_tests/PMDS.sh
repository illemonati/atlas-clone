#! /bin/bash

. $(dirname $0)/find_atlas

mkdir PMDS
file=PMDS/times
if [ ! -f "$file" ]; then
	echo -e "PMDS" > PMDS/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task PMDS --bam simulate/bamFile.bam --fasta simulate/bamFile.fasta --out PMDS/PMDS --logFile PMDS/PMDS.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> PMDS/times
