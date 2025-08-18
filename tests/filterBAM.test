#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --F 0.1 --fixedSeed 99 \
  --sampleSize 9 --chrLength 1111 --fracPoly 1.0 \
  --alpha 2.0 --beta 2.0 --seqType paired --seqCycles 101

echo "chr1	0	1000" > bed.bed
echo "chr2	100	800" >> bed.bed
echo "chr3	600	760" >> bed.bed

out="filterBAM_0"
$atlas --task filterBAM --dryRun \
           --bam simulate_ind1.bam --regions bed.bed \
           --fixedSeed 91 --out $out --logFile $out.out 2> $out.eout

out="filterBAM_1"
$atlas --task filterBAM --dryRun \
	   --bam simulate_ind1.bam --mask bed.bed --maskPorosity 0.1 \
	   --fixedSeed 92 --out $out --logFile $out.out 2> $out.eout

out="filterBAM_2"
$atlas --task filterBAM --makeSingle --bam simulate_ind1.bam \
	   --fixedSeed 93 --out $out --logFile $out.out 2> $out.eout

out="bd_2"
$atlas --task BAMDiagnostics --bam filterBAM_2_filtered.bam  \
	   --fixedSeed 94 --out $out --logFile $out.out 2> $out.eout

for i in {3..9}; do
	samtools view simulate_ind$i.bam | head -250 | tail -10 | cut -f1 \
		   > blacklist_$i.txt
	u=$(echo "$i*10" | bc)
	out="filterBAM_$i"
	$atlas --task filterBAM --bam simulate_ind$i.bam \
		   --filterMQ 0,$u --blacklist blacklist_$i.txt --filterReadLength 0,$u \
		   --filterFragmentLength 0,$u --filterSoftClips "0.$i" \
		   --fixedSeed 1$i --out $out --logFile $out.out 2> $out.eout
done

