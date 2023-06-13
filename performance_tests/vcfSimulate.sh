#! /bin/bash

. $(dirname $0)/find_atlas

mkdir simulate
file=simulate/VCFtimes
if [ ! -f "$file" ]; then
	echo -e "simulateVCF" > simulate/VCFtimes
fi

timeFor10Runs=0
for i in {1..5}; do
start=`date +%s.%N`

$atlas --task simulate --type HW --sampleSize 10 --chrLength 500000{2},300000,400000,600000 --ploidy 2{3},1,2 --depth 10,8{2},5{2} --baseFreq 0.5,0.3,0.2,0 --refDiv 0.5 --vcf --fixedSeed 0 --out simulate/vcfFile

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> simulate/VCFtimes

