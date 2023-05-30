#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 10 --chrLength 1000 --ploidy 2

$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format beagle --out ATLAS_beagle --logFile convertVCF_beagle.out
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format geno --out ATLAS_geno --logFile convertVCF_geno.out
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format LFMM --genotypes call --out ATLAS_LFMM_call --logFile convertVCF_LFMM_call.out
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format LFMM --genotypes posterior --maxMissing 0.01 --out ATLAS_LFMM_posterior --logFile convertVCF_LFMM_posterior.out
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format posfile --out ATLAS_posfile --logFile convertVCF_posfile.out
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format genfile --out ATLAS_genfile --logFile convertVCF_genfile.out
