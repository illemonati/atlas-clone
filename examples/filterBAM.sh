#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --F 0.1 --fixedSeed 99 \
  --sampleSize 19 --chrLength 1111 --fracPoly 1.0 \
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
	   --fixedSeed 91 --out $out --logFile $out.out 2> $out.eout

for i in {2..19}; do
	samtools view simulate_ind$i.bam | head -250 | tail -10 | cut -f1 \
		   > blacklist_$i.txt
	u=$(echo "$i*5" | bc)
	out="filterBAM_$i"
	$atlas --task filterBAM --bam simulate_ind$i.bam \
		   --filterMQ 0,$u --blacklist blacklist_$i.txt --filterReadLength 0,$u \
		   --filterFragmentLength 0,$u --filterSoftClips "0.$i" \
		   --fixedSeed 1$i --out $out --logFile $out.out 2> $out.eout
done

