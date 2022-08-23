#! /bin/bash

. $(dirname $0)/find_atlas

N=20
if [[ "$#" -eq 1 ]]; then
    N=$1
fi

echo "doing $N samples"

$atlas --task simulate --type HW --F 0.1 --chrLength 1000 --sampleSize $N --fracPoly 1.0 --alpha 2.0 --beta 2.0 --fixedSeed 0 --logFile simulate.out

for f in *.bam; do
    $atlas --task GLF --bam $f --fasta ATLAS_simulations.fasta --fixedSeed 0 --logFile GLF_${f%.bam}.out
done

allSamples=`find . -path '*_ind*.glf.gz' | paste -s -d ',' -`

$atlas --task majorMinor --glf $allSamples --fixedSeed 0 --logFile majorMinor.out
