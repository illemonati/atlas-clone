#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44 --chrLength 1000{3} --depth 0

######################################################
# Test Bed2Fastq
######################################################

# write bed file for testing
(
	echo "chr1 200 201 Info1"
	echo "chr1 400 401 Info2"
	echo "chr2 50 51 Info3"
	echo "chr2 950 951 Info4"
	echo "chr3 200 210 Info5"
) > liftOver.bed 

out="bed2fastq"
flank=123
delim="___-___"
$atlas --task liftOver --mode Bed2Fastq --bed liftOver.bed --fasta simulate.fasta --flank $flank\
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout
	   
# compare read names in fastq (chrLength = 1000, see above)
awk -vflank=$flank -vdelim=$delim '{for(i=$2; i<$3; ++i){ start=i-flank; if(start < 0){start=0}; end=i+flank+1; if(end>1000){end = 1000}; print "@" $1 ":" start "-" end delim i-start+1 delim $4}}' liftOver.bed > expectedNames.txt
awk '$1~/@/' ${out}.fastq > atlasNames.txt
diff expectedNames.txt atlasNames.txt >> $out.eout

# compare sequences
awk -vflank=$flank '{for(i=$2; i<$3; ++i){ start=i-flank+1; if(start < 1){start=1}; end=i+flank+1; if(end>1000){end = 1000}; print $1 ":" start "-" end}}' liftOver.bed | while read line
do
	samtools faidx -n 1000 simulate.fasta $line
done | awk '$1 !~/^>/' > expectedSequences.txt 
awk '{l++}$1~/^@/{l=1} l==2{print $0}' ${out}.fastq > atlasSequences.txt
diff expectedSequences.txt atlasSequences.txt >> $out.eout

######################################################
# Test Bam2Bed
######################################################

samtools view -H simulate.bam > testBam.sam # header
( 
	echo "ch1:195-206${delim}6${delim}Info1 0 chr1 195 50 11M * 0 11 AAAAAAAAAAA 12345678901"  
	echo "ch1:395-406${delim}6${delim}Info2 0 chr1 395 50 3S6M2S * 0 11 AAAAAAAAAAA 12345678901" # softclipping
	echo "ch1:400-411${delim}6${delim}notMapped 0 chr1 400 50 10S1M * 0 11 AAAAAAAAAAA 12345678901" # not mapped at pos bc softclipping
	echo "ch1:400-411${delim}6${delim}notMapped 0 chr1 400 50 5M1I5M * 0 11 AAAAAAAAAAA 12345678901" # not mapped at pos bc insertion
	echo "ch2:45-56${delim}6${delim}Info3 0 chr2 45 50 11M * 0 11 AAAAAAAAAAA 12345678901"  
	echo "ch2:945-956${delim}6${delim}Info4 0 chr2 945 50 11M * 0 11 AAAAAAAAAAA 12345678901"
	
) >> testBam.sam

for i in ${0..9}; do
	start=$(echo 195 + $i | bc)
	end=$(echo 206 + $i | bc)
	echo "ch3:${start}-${end}${delim}6${delim}Info5 0 chr3 $start 50 11M * 0 11 AAAAAAAAAAA 12345678901" >> testBam.sam
done

samtools view -bS testBam.sam > testBam.bam

$atlas --task liftOver --mode Bam2Bed --bam testBam.bam \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

# the initial bed file (liftOver.bed) should be the same as the atlas output ($out.bed)
diff  liftOver.bed ${out}.bed >> $out.eout
