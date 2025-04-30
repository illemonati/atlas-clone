#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44 --chrLength 1000{3} --depth 0
echo "$atlas"

######################################################
# Test Bed2Fastq
######################################################

# write bed file for testing
(
	echo "chr1 200 201 Info1"
	echo "chr1 400 401 softclipping"
	echo "chr1 420 421 mapping4M1I6M"
	echo "chr1 430 431 mapping6M1I4M"
	echo "chr2 50 51 softclippingLeft"
	echo "chr2 850 851 softclippingRight"
	echo "chr2 950 951 deletion"
	echo "chr3 200 210 Info2"
) > liftOver.bed 

out="bed2fastq"
flank=123
delim="___-___"

# running mode 1: Bed2Fastq
$atlas liftOver --mode Bed2Fastq --bed liftOver.bed --fasta simulate.fasta --flank $flank\
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout
	   
# compare read names in fastq (chrLength = 1000, see above)
awk -vflank=$flank -vdelim=$delim '{for(i=$2; i<$3; ++i){ start=i-flank; if(start < 0){start=0}; end=i+flank+1; if(end>1000){end = 1000}; print "@" $1 ":" start "-" end delim i-start delim $4}}' liftOver.bed > expectedNames.txt
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
	echo "ch1:195-206${delim}5${delim}Info1 0 chr1 196 50 11M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1"  
	echo "ch1:395-406${delim}5${delim}softclipping 0 chr1 396 50 3S6M2S * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # softclipping
	echo "ch1:405-416${delim}5${delim}notMappedSoftClipping 0 chr1 406 50 10S1M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # not mapped at pos bc softclipping
	echo "ch1:415-426${delim}5${delim}mapping4M1I6M 0 chr1 417 50 4M1I6M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # mapped with insertion
	echo "ch1:425-436${delim}5${delim}notMappedInsertion 0 chr1 426 50 5M1I5M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # not mapped at pos bc insertion
	echo "ch1:435-446${delim}5${delim}mapping6M1I4M 0 chr1 426 50 6M1I4M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # mapped with insertion
	echo "ch2:45-56${delim}5${delim}softclippingLeft 0 chr2 46 50 2S9M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1"  # softclipping at start
	echo "ch2:845-956${delim}5${delim}softclippingRight 0 chr2 846 50 8M3S * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # softclipping at end
	echo "ch2:955-966${delim}5${delim}deletion 0 chr2 945 50 5M1D6M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" # deletion
	
) | tr ' ' '\t' >> testBam.sam


for i in {0..9}; do
	start=$(echo 195 + $i | bc)
	end=$(echo 206 + $i | bc)
	echo "ch3:${start}-${end}${delim}6${delim}Info2 0 chr3 $start 50 11M * 0 11 AAAAAAAAAAA 12345678901 RG:Z:SimReadGroup1" | tr ' ' '\t' >> testBam.sam
done

samtools view -bS testBam.sam > testBam.bam

# running mode 2: Bam2Bed
$atlas liftOver --mode Bam2Bed --bam testBam.bam \
	   --fixedSeed 42 --out $out --logFile $out.out 2>> $out.eout


# the initial bed file (liftOver.bed) after some formatting should be the same as the atlas output ($out.bed)
awk 'BEGIN {OFS="\t"} {if ($2+1==$3) {print $1, $2, $3, $4} else {pos0=$2; while (pos0 < $3) {print $1, pos0, pos0+1, $4; pos0++}}}' liftOver.bed > expectedBed.bed
diff expectedBed.bed ${out}.bed >> $out.eout
