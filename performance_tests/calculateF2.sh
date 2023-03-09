#! /bin/bash

. $(dirname $0)/find_atlas

mkdir calculateF2
file=calculateF2/times
if [ ! -f "$file" ]; then
	echo -e "calculateF2" > calculateF2/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task calculateF2 --vcf simulate/vcfFile.vcf.gz --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 --logFile calculateF2/calculateF2.out --out calculateF2/out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> calculateF2/times
