#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}

out=$name/${name}_beagle
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --chr chr1 --format beagle --out $out --logFile $out.out 2> $out.err > /dev/null

out=$name/${name}_geno
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --format geno --out $out --logFile $out.out 2> $out.err > /dev/null

out=$name/${name}_call
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --format LFMM --genotypes call --out $out --logFile $out.out 2> $out.err > /dev/null

out=$name/${name}_post
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --format LFMM --genotypes posterior --maxMissing 0.01 --out $out --logFile $out.out 2> $out.err > /dev/null

out=$name/${name}_posfile
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --format posfile --out $out --logFile $out.out 2> $out.err > /dev/null

out=$name/${name}_genfile
$atlas --task convertVCF --vcf VCF/VCF.vcf.gz --format genfile --out $out --logFile $out.out 2> $out.err > /dev/null
