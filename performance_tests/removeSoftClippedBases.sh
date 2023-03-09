#! /bin/bash

. $(dirname $0)/find_atlas

mkdir removeSoftClippedBases
file=removeSoftClippedBases/times
if [ ! -f "$file" ]; then
	echo -e "removeSoftClippedBases" > removeSoftClippedBases/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task removeSoftClippedBases --bam simulate/bamFile.bam --fixedSeed 0 --out removeSoftClippedBases/removeSoftClippedBases --logFile removeSoftClippedBases/removeSoftClippedBases.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> removeSoftClippedBases/times
