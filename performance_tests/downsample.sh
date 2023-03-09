#! /bin/bash

. $(dirname $0)/find_atlas

mkdir downsample
file=downsample/times
if [ ! -f "$file" ]; then
	echo -e "downsample" > downsample/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task downsample --bam simulate/bamFile.bam --prob 0.5,0.3,0.1 --fixedSeed 0 --logFile downsample/downsample.out --out downsample/downsampled


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> downsample/times
