#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

N=20
if [[ "$#" -eq 1 ]]; then
    N=$1
fi

echo "doing $N samples"

$atlas --task simulate --type BAM:HW --F 0.1 --chrLength 1000 --sampleSize $N --fracPoly 1.0 --alpha 2.0 --beta 2.0 --fixedSeed 0

for f in *.bam; do
    $atlas --task GLF --bam $f --fasta ATLAS_simulations.fasta --fixedSeed 0
done

allSamples=`find . -path '*_ind*.glf.gz' | paste -s -d ',' -`

$atlas --task majorMinor --glf $allSamples --fixedSeed 0
