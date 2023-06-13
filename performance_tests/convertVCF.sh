#! /bin/bash

. $(dirname $0)/find_atlas

mkdir convertVCF
file=convertVCF/times
if [ ! -f "$file" ]; then
	echo -e "convertVCF" > convertVCF/times
fi

timeFor10Runs=0

for i in {1..2}; do
start=`date +%s.%N`

$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format beagle --out convertVCF/ATLAS_beagle --logFile convertVCF/convertVCF_beagle.out
$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format geno --out convertVCF/ATLAS_geno --logFile convertVCF/convertVCF_geno.out
$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format LFMM --genotypes call --out convertVCF/ATLAS_LFMM_call --logFile convertVCF/convertVCF_LFMM_call.out
$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format LFMM --genotypes posterior --maxMissing 0.01 --out convertVCF/ATLAS_LFMM_posterior --logFile convertVCF/convertVCF_LFMM_posterior.out
$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format posfile --out convertVCF/ATLAS_posfile --logFile convertVCF/convertVCF_posfile.out
$atlas --task convertVCF --vcf simulate/vcfFile.vcf.gz --format genfile --out convertVCF/ATLAS_genfile --logFile convertVCF/convertVCF_genfile.out


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> convertVCF/times
