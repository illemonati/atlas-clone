#include "TPositionBasedLiftOver.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/fillContainer.h"
#include "genometools/TFastaReader.h"
#include "genometools/TBed.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"

#include <algorithm>

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace PopulationTools {

namespace impl{
    const std::string Bed2FastqMode = "Bed2Fastq";
    const std::string Bam2BedMode = "Bam2Bed";
    const std::string readNameDelimiter = "___-___";
}
    

//--------------------------------------------
// TBedToFastq
//--------------------------------------------

void TBedToFastq::run(){      
    // reading user input
    const auto fastaName = parameters().get("fasta"); 
    logfile().listFlush("Will read reference from FASTA file '", fastaName, "'. (parameter 'fasta')");
    genometools::TFastaReader reference(fastaName); 
    const genometools::TChromosomes& chrs = reference.chromosomes();
    
    const auto bedname = parameters().get("bed");
    logfile().startIndent("Reading windows to parse (parameter 'bed'):");
    genometools::TBedWithInfo<std::string> bed(bedname, chrs); // reading correct refID from reference
    logfile().endIndent();

    const long flank = parameters().get("flank", 100);
    logfile().write("Will write reads extendeding ", flank, " bases on both sides. (paramter 'flank').");
    const auto outname = parameters().get("out", "ATLAS_" + impl::Bed2FastqMode);
    logfile().write("Will write output FASTQ to file '", outname, "'. (parameter 'out')");

    // open output file
    coretools::TOutputFile out(outname + ".fastq", "");
    
    coretools::TTimer timer;
    constexpr size_t dCounter = 100000;
    size_t counter            = 0;    
    size_t nextPrint          = dCounter;

    // looping over windows
    logfile().startIndent("Writing reads to FASTQ file:");
    for (const auto& b: bed) {
        // for every position within a range, write one fastq sequence

        for (size_t i = 0; i < b.size(); ++i) {                     
            const long start = std::max((long) 0, (long) b.fromOnChr() + (long) i - flank); // included
            const size_t end = std::min((size_t) chrs[b.refID()].length(), b.fromOnChr() + i + flank + 1); // not included
            const size_t posInRead = b.fromOnChr() + i - start + 1;
            const size_t seqLength = end - start;

            // write
            out.writeln("@", chrs[b.refID()].name(), ":", start, "-", end, impl::readNameDelimiter, posInRead, impl::readNameDelimiter, b.info()); // first thing we actually need for the fastq, but we need to add more info
            out.writeln(reference.view(b.refID(), start, seqLength));
            out.writeln("+");
            out.writeln(std::string(seqLength, 'F'));
        }
        
        // report progress
        if (counter >= nextPrint) {
            logfile().list("Parsed ", counter, " positions in ", timer.formattedTime(), ".");
            nextPrint += dCounter;
        }
    }
    logfile().list("Parsed ", counter, " positions in ", timer.formattedTime(), ".");    
}


//--------------------------------------------
// TBamToBed
//--------------------------------------------

void TBamToBed::_handleAlignment(BAM::TAlignment& alignment){
    

    // read position to extrac from name    
    std::vector<std::string_view> vec;
    coretools::str::fillContainerFromString(alignment.name(), vec, impl::readNameDelimiter);

    if(vec.size() != 3){
        UERROR("Unable to parse name of alignment '", alignment.name(), "': did you map a FASTQ file produced with task=BLAH mode=", impl::Bed2FastqMode, "?");
    }
    
    size_t posInRead = coretools::str::fromString<size_t>(vec[1]);
    if(alignment.isAlignedAtInternalPos(posInRead)){
        _outBed.add(alignment.positionInRef(posInRead), vec[2]);
    } else {
        ++_numPosNotAligned;
    }    
}

void TBamToBed::run(){
    logfile().list("Extracting BED information from mapped reads (task=BLAH, mode=", impl::Bam2BedMode, ")");

    // open output BED    
    std::string outname = parameters().get("out", "ATLAS_" + impl::Bam2BedMode);
    outname += ".bed";
    logfile().write("Will write output BED to file '", outname, "'. (parameter 'out')");

    _outBed.setChromosomes(_genome.bamFile().chromosomes());
    _numPosNotAligned = 0;

	// traverse BAM
	_traverseBAMPassedQC();

    // report and write BED    
    _outBed.write(outname);
    logfile().conclude("Wrote a total of ", _outBed.size(), " positions.");
    if(_numPosNotAligned > 0){
        logfile().conclude("Ignored ", _numPosNotAligned, ", positions that did not align.");
    }
}

//--------------------------------------------
// TPositionBasedLiftOver
//--------------------------------------------

void TPositionBasedLiftOver::run(){
    const auto mode = parameters().get("mode");
    if(mode == impl::Bed2FastqMode){
        TBedToFastq converter;
        converter.run();
    } else if (mode == impl::Bam2BedMode){

    } else {
        UERROR("Unknown mode '", mode, "'! Accepted values are '", impl::Bed2FastqMode, "', and '", impl::Bam2BedMode, "'.");
    }
}

} // end namespace poptools
