#! /bin/bash

. $(dirname $0)/find_atlas

mkdir readOverlap
file=readOverlap/times
if [ ! -f "$file" ]; then
	echo -e "readOverlap" > readOverlap/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task readOverlap --bam simulate/bamFile.bam --fixedSeed 0 --out readOverlap/readOverlap --logFile readOverlap/readOverlap.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> readOverlap/times
