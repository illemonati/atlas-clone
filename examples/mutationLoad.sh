#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out --chrLength 100000

# Create a bed file with every 7th position
printf "chr1\t0\t1\n" > bed.bed
for i in {7..10000..7}; do
	printf "chr1\t$i\t$((i+1))\n" >> bed.bed
done

# Calculate mutation load
$atlas mutationLoad --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --bed bed.bed --logFile mutationLoad.out
