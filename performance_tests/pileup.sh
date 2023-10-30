#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task pileup --readUpToDepth 100 --printAll \
	   --bam BAM/BAM.bam --fields depth,alleles,mates,strands \
	   --window 100000 --histograms depth,allelicDepth,contexts,qualities \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null
