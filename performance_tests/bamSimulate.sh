#! /bin/bash

. $(dirname $0)/find_atlas

mkdir simulate
file=simulate/BAMtimes
if [ ! -f "$file" ]; then
	echo -e "simulateBAM" > simulate/BAMtimes
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task simulate --chrLength 500000{2},300000,400000,600000 --ploidy 2{3},1,2 --depth 10,8{2},5{2} --baseFreq 0.5,0.3,0.2,0 --refDiv 0.5 --seqType paired --seqCycles 200 --fragmentLength 'fixed(500)' --baseQuality 'binomial(95,0.01)[0,93]' --mappingQuality 'normal(60,10)[1,255]' --softClipping 'normal(10,10)[0,20]' --pmd 'doubleStrand:Exponential[10,20,1,2]:Empiric[0.3]' --recal "intercept[0.1];quality:polynomial[0.8,0.2,0.01]" --frequency 0.2 --fixedSeed 0 --out simulate/bamFile

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> simulate/BAMtimes

