#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type SFS --theta 0.001 --fixedSeed 212 \
  --chrLength 15111 --depth 11 --ploidy 2 --sampleSize 1

name="GLF"
samples=""
for f in *.bam; do
	out=${f%.bam}
	$atlas --task GLF --bam $f  \
	   --fixedSeed 215 --out $out --logFile $out.out 2> $out.eout
	samples="${samples}$out.glf.gz,"
done

out="saf"
$atlas --task saf --glf $samples --fasta simulate.fasta \
	   --fixedSeed 219 --out $out --logFile $out.out 2> $out.eout
