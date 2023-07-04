#!/bin/bash

#simulate both Fastq and BAM files
#at least two individuals
sampleSize=3
chrLength=5000

mkdir Simulations
cd Simulations
../../build/atlas simulate --fastq,bam --type HW --sampleSize $sampleSize --chrLength $chrLength --fracPoly 0.05 --alpha 0.3 --beta 0.8 --depth 5


# Set the source directory path
source_directory="./"

# Set the destination directory path
destination_directory="./"

# Array of file extensions to check
extensions=("bam" "fastq")

# Iterate over the files in the source directory
for file in "$source_directory"/*; do
    # Check if the item is a file
    if [[ -f "$file" ]]; then
        # Get the file extension
        extension="${file##*.}"
        
        # Check if the extension is in the list
        if [[ " ${extensions[*]} " == *" $extension "* ]]; then
            # Create the destination subdirectory if it doesn't exist
            destination_subdirectory="$destination_directory/$extension"
            if [[ ! -d "$destination_subdirectory" ]]; then
                mkdir -p "$destination_subdirectory"
            fi
            
            # Move the file to the destination subdirectory
            mv "$file" "$destination_subdirectory/$(basename "$file")"
            #echo "Moved file from: $(basename "$file") to $extension/"
        fi
    fi
done

echo
echo "#########################################################"
echo "Files simulated and sorted in respective directories"
echo "#########################################################"
echo

#check fastq files against bam files converted to fastq (this direction of comparison erases the problem of different read groups in bam file that would be seen as missmatches or non existing in previous fastq file)

mkdir catFastq
cd ./fastq

#create fastq files not sorted per RG -> _0/_1 should not exist anymore
#create fastq according to the various individuals name

for ((i=1;i<=sampleSize;i++)); do
	filename="ATLAS_simulations_ind${i}";
	touch $filename.fastq ;
        #echo "${filename}.fastq" ;
        #echo
    #for every new individual, concatenate its RG into this newly created file
    cat ATLAS_simulations_ind${i}_* >>  ATLAS_simulations_ind${i}.fastq 
    mv ATLAS_simulations_ind${i}.fastq ../catFastq
done

echo
cd ../catFastq
echo

#checks the content of every catfastq file against the bam2fastq (using grep) in order to see how many lines do not match
#additionally create a percentage error rate to give a value to how many data is missing from catFastq to bam2fastq

sPattern="^[ACGT]+$"
mPattern="^@"
touch ../../missedSequencesOutput.txt
touch ../grepErrors.txt

for filex in *; do 
    fastqFile="$filex"
    bamFile="../bam/${filex%.fastq}.bam"
    
    totLines=0
    notSequenceLines=0
    notMetaLines=0
    notQualitiesLines=0

    echo
    echo "Checking lines for: $filex"
    while IFS= read -r line; do
    #if [[ "$line" =~ $pattern ]] ; then
    #    echo
    #    echo "$line"
    #    echo
    #fi

       if ! samtools view "$bamFile" | grep -Fqe "$line" >> ../grepErrors.txt 2>&1; then
            #echo "$line"
            if [[ "$line" =~ $sPattern ]] ; then
                notSequenceLines=$((notSequenceLines+1))
                #echo -e "$line \tNOT MATCHED"
            elif [[ "$line" =~ $mPattern ]] ; then
                notMetaLines=$((notMetaLines+1))
                #echo -e "METADATA \tNOT MATCHED"
            else
                notQualitiesLines=$((notQualitiesLines+1))
                echo "$line"
            fi
        fi
        totLines=$((totLines+1))
    done < "$fastqFile"

    #echo
    #echo -e "For ${fastqFile%.fastq}:
    #        \t the sequence(s) not found from fastq to bam are: $notSequenceLines ($((notSequenceLines*100/totLines))%);
    #        \t the metadata not found are: $notMetaLines ($((notMetaLines*100/totLines))%);
    #        \t the qualities not found are: $notQualitiesLines ($((notQualitiesLines*100/totLines))%);
    #        \t over a total of $totLines \n\n" >> ../missedSequencesOutput.txt
    #echo "Finished writing $fastqFile results"
    #echo
    quality=$((1-(notSequenceLines+notQualitiesLines)/totLines))
    echo -e "For ${fastqFile%.fastq}:
            \t the sequence(s) not found from fastq to bam are: $notSequenceLines ($((notSequenceLines*100/totLines))%);
            \t the qualities not found are: $notQualitiesLines ($((notQualitiesLines*100/totLines))%);
            \t over a total of $totLines
            \t this means a quality of transcription of: ($((quality * 100)))% \n\n" >> ../../missedSequencesOutput.txt
    echo "Finished writing $fastqFile results"
    echo
done
