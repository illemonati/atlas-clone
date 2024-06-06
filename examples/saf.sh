#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type SFS --theta 0.001 --chrLength 15111 --depth 11 --ploidy 2 --sampleSize 1

name="GLF"
samples=""
for f in *.bam; do
	name=${f%.bam}
	$atlas --task GLF --bam $f  \
	   --fixedSeed 0 --out $name --logFile glf_$name.out
	samples="${samples}$name.glf.gz,"
done

name="saf"
$atlas --task saf --glf $samples --fasta ATLAS_simulations.fasta \
	   --fixedSeed 0 --out $name --logFile $name.out
