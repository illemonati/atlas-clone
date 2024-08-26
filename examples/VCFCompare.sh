#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas

. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 2321 --ploidy 2 \
	--out simulate1 --logFile simulate1.out --fixedSeed 1

. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 2321 --ploidy 2 \
	--out simulate2 --logFile simulate2.out --fixedSeed 2

out="VCFCompare_f_ss"
$atlas --task VCFCompare --vcf simulate1.vcf.gz --samples ind_0,ind_1 \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

out="VCFCompare_ff_s"
$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz --samples ind_0 \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

out="VCFCompare_ff_ss"
$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz \
	   --samples ind_0,ind_1 \
	   --fixedSeed 5 --out $out --logFile $out.out 2> $out.eout
