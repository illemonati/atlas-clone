#! /bin/bash

. $(dirname $0)/find_atlas

. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 2321 --ploidy 2 \
	--out simulate1 --logFile simulate1.out --fixedSeed 1

. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 2321 --ploidy 2 \
	--out simulate2 --logFile simulate2.out --fixedSeed 2

$atlas --task VCFCompare --vcf simulate1.vcf.gz --samples ind_0,ind_1 \
	   --fixedSeed 0 --out VCFCompare_f_ss --logFile VCFCompare_f_ss.out

$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz --samples ind_0 \
	   --fixedSeed 0 --out VCFCompare_ff_s --logFile VCFCompare_ff_s.out

$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz --samples ind_0,ind_1 \
	   --fixedSeed 0 --out VCFCompare_ff_ss --logFile VCFCompare_ff_ss.out
