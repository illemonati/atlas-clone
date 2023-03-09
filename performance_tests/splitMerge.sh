#! /bin/bash

. $(dirname $0)/find_atlas

mkdir splitMerge
file=splitMerge/times
if [ ! -f "$file" ]; then
	echo -e "splitMerge" > splitMerge/times
fi

timeFor10Runs=0
echo "readGroup seqType seqCycles" > splitMerge/rgs_paired.txt
echo "SimReadGroup1 paired 200" >> splitMerge/rgs_paired.txt

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task splitMerge --bam simulate/bamFile.bam --readGroupSettings splitMerge/rgs_paired.txt --logFile splitMerge/splitMerge.out --mergingMethod middle --fixedSeed 0 --out splitMerge/splitMerge

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> splitMerge/times
