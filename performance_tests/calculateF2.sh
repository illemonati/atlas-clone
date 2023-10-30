#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task calculateF2 \
	   --vcf VCF/VCF.vcf.gz --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 \
	   --out $out --logFile $out.out 2> $out.err > /dev/null
