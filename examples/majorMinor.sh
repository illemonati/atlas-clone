#! /bin/bash

. $(dirname $0)/find_atlas

N=97
if [[ "$#" -eq 1 ]]; then
    N=$1
fi

echo "doing $N samples"

. $(dirname $0)/simulate --type HW --F 0.1 --chrLength 1432 \
	--sampleSize $N --fracPoly 0.1 --alpha 2.0 --beta 2.0

for f in *.bam; do
	name=${f%.bam}
    $atlas --task GLF --bam $f --fasta ATLAS_simulations.fasta \
		   --fixedSeed 0 --out GLF_$name --logFile GLF_$name.out
done

allSamples=`find . -path '*_ind*.glf.gz' | paste -s -d ',' -`

name="majorMinor"
$atlas --task majorMinor --glf $allSamples --method Skotte \
	   --minMAF 0.05 --maxThreads 1 --bgz --minSamplesWithData 83 \
	   --fixedSeed 0 --out $name --logFile $name.out

name="Skotte_fasta"
$atlas --task majorMinor --method Skotte --minSamplesWithData 83 \
	--glf $allSamples --fasta ATLAS_simulations.fasta \
	--minMAF 0.05 --maxThreads 1 --bgz \
	--fixedSeed 0 --out $name --logFile $name.out

name="MLE_fasta"
$atlas --task majorMinor --method MLE --minSamplesWithData 83 \
	--glf $allSamples --fasta ATLAS_simulations.fasta \
	--minMAF 0.05 --maxThreads 1 --bgz \
	--fixedSeed 0 --out $name --logFile $name.out

echo "chr pos ref alt" > alleles.txt
zcat Skotte_fasta.vcf.gz | tail -n +10 | awk '{if (NR % 3 == 0) {print $1, $2, $4, $5}}' >> alleles.txt

name="alleles"
$atlas --task majorMinor --method Skotte --minSamplesWithData 83 \
	--glf $allSamples --alleles alleles.txt \
	--minMAF 0.05 --maxThreads 1 \
	--fixedSeed 0 --out $name --logFile $name.out
