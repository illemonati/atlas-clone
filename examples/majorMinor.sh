#! /bin/bash

. $(dirname $0)/find_atlas

N=97
if [[ "$#" -eq 1 ]]; then
    N=$1
fi

echo "doing $N samples"

. $(dirname $0)/simulate --type HW --F 0.1 --chrLength 1432 \
	--sampleSize $N --fracPoly 0.1 --alpha 2.0 --beta 2.0 --fixedSeed 133

for f in *.bam; do
	out=GLF_${f%.bam}
    $atlas --task GLF --bam $f --fasta simulate.fasta \
		   --fixedSeed 131 --out $out --logFile $out.out 2> $out.eout
done

allSamples=`find . -path '*_ind*.glf.gz' | paste -s -d ',' -`

out="majorMinor"
$atlas --task majorMinor --glf $allSamples --method Skotte \
	   --minMAF 0.05 --maxThreads 1 --bgz --minSamplesWithData 83 \
	   --fixedSeed 132 --out $out --logFile $out.out 2> $out.eout

out="Skotte_fasta"
$atlas --task majorMinor --method Skotte --minSamplesWithData 83 \
	   --glf $allSamples --fasta simulate.fasta \
	   --minMAF 0.05 --maxThreads 1 --bgz \
	   --fixedSeed 134 --out $out --logFile $out.out 2> $out.eout

out="MLE_fasta"
$atlas --task majorMinor --method MLE --minSamplesWithData 83 \
	   --glf $allSamples --fasta simulate.fasta \
	   --minMAF 0.05 --maxThreads 1 --bgz \
	   --fixedSeed 135 --out $out --logFile $out.out 2> $out.eout

echo "chr pos ref alt" > alleles.txt
gunzip -c Skotte_fasta.vcf.gz | tail -n +10 | awk '{if (NR % 3 == 0) {print $1, $2, $4, $5}}' >> alleles.txt

out="alleles"
$atlas --task majorMinor --method Skotte --minSamplesWithData 83 \
	   --glf $allSamples --alleles alleles.txt \
	   --minMAF 0.05 --maxThreads 1 \
	   --fixedSeed 136 --out $out --logFile $out.out 2> $out.eout
