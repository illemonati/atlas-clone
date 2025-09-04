#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate first BAM File
$atlas simulate --logFile first.out --out first

# Simulate three BAM Files using fasta-file of first bam file
$atlas simulate --fasta first.fasta --sampleSize 3 --logFile multi.out --out multi
