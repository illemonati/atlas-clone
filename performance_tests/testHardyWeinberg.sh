#! /bin/bash

. $(dirname $0)/find_atlas

mkdir testHardyWeinberg
file=testHardyWeinberg/times
if [ ! -f "$file" ]; then
	echo -e "testHardyWeinberg" > testHardyWeinberg/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task testHardyWeinberg --vcf simulate/vcfFile.vcf.gz --out testHardyWeinberg/testHardyWeinberg --logFile testHardyWeinberg/testHardyWeinberg.out 


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> testHardyWeinberg/times
