#! /bin/bash

. $(dirname $0)/find_atlas

mkdir filterBAM
file=filterBAM/times
if [ ! -f "$file" ]; then
	echo -e "filterBAM" > filterBAM/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

u=$(echo "$i*10" | bc)
$atlas --task filterBAM --bam simulate/bamFile.bam --filterSoftClips --fixedSeed 0 --filterMQ 0,$u --filterReadLength 0,$u --filterFragmentLength 0,$u --logFile filterBAM/filterBAM.out


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> filterBAM/times
