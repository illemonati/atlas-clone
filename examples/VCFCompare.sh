#! /bin/bash

. $(dirname $0)/find_atlas


$atlas --task simulate --vcf --type HW --sampleSize 2 --chrLength 2000 --ploidy 2 --fixedSeed 1 --out simulate1 --logFile simulate1.out
$atlas --task simulate --vcf --type HW --sampleSize 2 --chrLength 2000 --ploidy 2 --fixedSeed 2 --out simulate2 --logFile simulate2.out

$atlas --task VCFCompare --vcf simulate1.vcf.gz --samples ind_0,ind_1 --out VCFCompare_f_ss --logFile VCFCompare_f_ss.out
$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz --samples ind_0 --out VCFCompare_ff_s --logFile VCFCompare_ff_s.out
$atlas --task VCFCompare --vcf simulate1.vcf.gz,simulate2.vcf.gz --samples ind_0,ind_1 --out VCFCompare_ff_ss --logFile VCFCompare_ff_ss.out
