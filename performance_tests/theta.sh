#! /bin/bash

. $(dirname $0)/find_atlas

mkdir theta
file=theta/times
if [ ! -f "$file" ]; then
	echo -e "theta" > theta/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

# theta
$atlas --task theta --bam simulate/bamFile.bam --printAll --fixedSeed 0 --extraVerbose --out theta/standard --logFile theta/standard.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> theta/times
