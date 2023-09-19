#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --writeTrueGenotypes

rm ATLAS_simulations.bam.bai # will automatically recreate it

# prepare alleles file
zcat ATLAS_simulations_trueGenotypes.vcf.gz | awk 'BEGIN{print "Chr", "Pos", "Allele1"}$1!~/#/{print $1, $2, $4}' > alleles.txt

# estimate mutation load
$atlas --task mutationLoad --bam ATLAS_simulations.bam --alleles alleles.txt --out mutationLoad --logFile mutationLoad.out 
