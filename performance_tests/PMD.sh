#! /bin/bash

. $(dirname $0)/find_atlas

mkdir PMD
file=PMD/times
if [ ! -f "$file" ]; then
	echo -e "PMD" > PMD/times
fi

timeFor10Runs=0

for i in {1..10}; do
start=`date +%s.%N`

$atlas --task PMD --bam simulate/bamFile.bam --fasta simulate/bamFile.fasta --pmdModels "doubleStrand:Exponential:Empiric"  --fixedSeed 0 --out PMD/PMD --logFile PMD/PMD.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> PMD/times
