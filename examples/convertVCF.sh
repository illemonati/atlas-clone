#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --vcf --type HW --sampleSize 10 --chrLength 100000 --fixedSeed 0

$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format beagle
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format geno
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format LFMM --genotypes call
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format LFMM --genotypes posterior --minSamplesWithData 10
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format posfile
$atlas --task convertVCF --vcf ATLAS_simulations.vcf.gz --format genfile
