#! /bin/bash

. $(dirname $0)/find_atlas

mkdir assessSoftClipping
file=assessSoftClipping/times
if [ ! -f "$file" ]; then
	echo -e "assessSoftClipping" > assessSoftClipping/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task assessSoftClipping --bam simulate/bamFile.bam --fixedSeed 0 --out assessSoftClipping/assessSoftClipping --logFile assessSoftClipping/assessSoftClipping.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> assessSoftClipping/times
