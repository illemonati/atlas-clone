#! /bin/bash

. $(dirname $0)/find_atlas

mkdir pileup
file=pileup/times
if [ ! -f "$file" ]; then
	echo -e "pileup" > pileup/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task pileup --bam simulate/bamFile.bam --fields depth,alleles,mates,strands --window 100000 --histograms depth,allelicDepth,contexts,qualities --readUpToDepth 100 --printAll --logFile pileup/pileup.out --out pileup/out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> pileup/times
