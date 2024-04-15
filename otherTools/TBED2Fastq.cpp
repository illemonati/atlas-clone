
#include "TBED2Fastq.h"

namespace Other{

//--------------------------------------------
// TBED2Fastq
//--------------------------------------------

TBED2Fastq::TBED2Fastq(std::string BedFileName, std::string FastaFileName, std::string OutFileName, int FlankingLength){
    // open bed file
    using coretools::instances::logfile;
    logfile().listFlush("Reading BED file '", BedFileName, "' ...");
    _bed.add(BedFileName);
    logfile().done();

    // open FASTA
    logfile().listFlush("Reading reference from FASTA file '", FastaFileName, "' ...");
    _fasta.open(FastaFileName);
    logfile().done();

    // opening output fastq file
    logfile().listFlush("Opening output fastq file '", OutFileName, "' ...");
    _out.open(OutFileName);
    logfile().done();
}


void TBED2Fastq::run(){
    using coretools::instances::logfile;
    logfile().startIndent("Parsing through BED file:");

    for(auto& pos: _bed){
        logfile().list(pos.from());
    }
}



} // end namesapce