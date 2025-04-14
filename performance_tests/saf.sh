#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task saf --glf GLF/GLF.glf.gz --fasta VCF/VCF.fasta \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null
