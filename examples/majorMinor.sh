#! /bin/bash

. $(dirname $0)/find_atlas

N=97
if [[ "$#" -eq 1 ]]; then
    N=$1
fi

echo "doing $N samples"

. $(dirname $0)/simulate --type HW --F 0.1 --chrLength 1432 \
	--sampleSize $N --fracPoly 0.1 --alpha 2.0 --beta 2.0

for f in *.bam; do
	name=${f%.bam}
    $atlas --task GLF --bam $f --fasta ATLAS_simulations.fasta \
		   --fixedSeed 0 --out GLF_$name --logFile GLF_$name.out
done

allSamples=`find . -path '*_ind*.glf.gz' | paste -s -d ',' -`

$atlas --task majorMinor --glf $allSamples --method Skotte \
	   --minMAF 0.05 --maxThreads 1 --bgz \
	   --fixedSeed 0 --out majorMinor --logFile majorMinor.out 
