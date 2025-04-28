#include "TFastaToFastq.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TFastaReader.h"
#include "genometools/TBed.h"
#include "coretools/Files/TOutputFile.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "TBamFile.h"



using coretools::instances::logfile;
using coretools::instances::parameters;


namespace PopulationTools {

    void writeFastqEntry(coretools::TOutputFile& out, const genometools::impl::TGenomeWindowWithInfo<double>& b, const genometools::TFastaReader& reference, size_t start, size_t seq_length) {
    
        out.writeln("@", reference.chrName(b.refID()), "_", start, "_", start+seq_length, "_", b.info()); // first thing we actually need for the fastq, but we need to add more info
        out.writeln(reference.view(b.refID(), start, seq_length));
        out.writeln("+");
        out.writeln(std::string(seq_length, 'F'));
    }


    TReadCenterAlignedBase::TReadCenterAlignedBase() : TBamReadTraverser<ReadType::Parsed>() {};

    void TReadCenterAlignedBase::_handleAlignment(BAM::TAlignment& alignment) {

};
    void TFastaToFastq::run(){

        const auto mode=parameters().get("mode");
        if (mode != "toFastq" && mode != "toBed") {
            UERROR("Unknown mode '", mode, "'! Accepted values are 'toFastq' and 'toBed'.");
        }

        const auto fastaname=parameters().get("fasta"); 
        const auto flank=parameters().get("flank", 100);
        const auto outname=parameters().get("out", "ATLAS_fastaToFastq");

        genometools::TFastaReader reference(fastaname); 

        if (mode == "toFastq") {
        
            const auto bedname=parameters().get("bed"); 
            genometools::TBedWithInfo<double> bed(bedname, reference.chromosomes()); // reading correct refID from reference

            coretools::TOutputFile out(outname + ".fastq", "");


            coretools::TTimer timer;
            constexpr size_t dCounter = 100000;
            size_t counter            = 0;
            size_t skip_counter       = 0;
            size_t nextPrint          = dCounter;

            const size_t seq_length = 2*flank + 1;


            for (const auto& b: bed) {

                // for every position within a range, write one fastq sequence
                for (size_t i = 0; i < b.size(); ++i) {
                    const int start = b.fromOnChr() + i - flank;
                    if (start >= 0) {
                        writeFastqEntry(out, b, reference, start, seq_length); 
                    } else {
                        ++skip_counter;
                        continue; // instead of skipping, go for max bt 0 and start-flank? 

                    }
                }
                
                // report progress
                if (counter >= nextPrint) {
                    logfile().list("Parsed ", nextPrint, " positions in ", timer.formattedTime(), ".");
                    nextPrint += dCounter;
                }

            }

            logfile().list("Removed ", skip_counter, " position(s) at the beginning of a chromosome.");

        } else {
            
            const auto bamname=parameters().get("bam");
            genometools::TBedWithInfo<double> out;




            bam = std::make_unique<BAM::TBamFile>(bamname, 0); // not sure what the ID is
            BAM::TAlignment alignment;




            logfile().list("Parsed length: ", alignment.parsedLength()); // I want to check this at some point bc I'm not sure this will always be 201



            out.write(outname + ".bed");

        }


    }
}

