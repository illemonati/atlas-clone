#! /bin/bash

. $(dirname $0)/find_atlas

mkdir inbreeding
file=inbreeding/times
if [ ! -f "$file" ]; then
	echo -e "inbreeding" > inbreeding/times
fi

timeFor10Runs=0
for i in {1..3}; do
start=`date +%s.%N`

$atlas --task inbreeding --vcf simulate/vcfFile.vcf.gz --numBurnin 1 --iterations 100 --fixedSeed 0 --logFile inbreeding/inbreeding.out --out inbreeding/out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> inbreeding/times
