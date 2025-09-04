#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 11 --chrLength 1111 \
  --ploidy 2 --fixedSeed 45

out="beagle"
$atlas --task convertVCF --vcf simulate.vcf.gz --format beagle \
	   --fixedSeed 46 --out $out --logFile $out.out 2> $out.eout

out="geno"
$atlas --task convertVCF --vcf simulate.vcf.gz --format geno \
	   --fixedSeed 47 --out $out --logFile $out.out 2> $out.eout

out="LFMM_call"
$atlas --task convertVCF --vcf simulate.vcf.gz --format LFMM \
	   --genotypes call \
	   --fixedSeed 48 --out $out --logFile $out.out 2> $out.eout

out="LFMM_post"
$atlas --task convertVCF --vcf simulate.vcf.gz --format LFMM \
	   --genotypes posterior --maxMissing 0.01 \
	   --fixedSeed 49 --out $out --logFile $out.out 2> $out.eout

out="posfile"
$atlas --task convertVCF --vcf simulate.vcf.gz --format posfile \
	   --fixedSeed 50 --out $out --logFile $out.out 2> $out.eout

out="genofile"
$atlas --task convertVCF --vcf simulate.vcf.gz --format genfile \
	   --fixedSeed 51 --out $out --logFile $out.out 2> $out.eout
