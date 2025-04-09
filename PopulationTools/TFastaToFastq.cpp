#include "TFastaToFastq.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TFastaReader.h"
#include "genometools/TBed.h"
#include "coretools/Files/TOutputFile.h"
#include "genometools/GenomePositions/TGenomeWindow.h"
#include "genometools/GenomePositions/TGenomePosition.h"



using coretools::instances::logfile;
using coretools::instances::parameters;



namespace PopulationTools {
    void TFastaToFastq::run(){
        logfile().list("hello");
        const auto bedfile=parameters().get("bed"); // andreas
        const auto fastafile=parameters().get("fasta"); // andreas
        const size_t flank=parameters().get("flanking_length", 3);
        const std::string outname=parameters().get<std::string>("out", "ATLAS_fastaToFastq");


        genometools::TFastaReader reference(fastafile); // andreas
        genometools::TBed bed(bedfile); // andreas

        coretools::TOutputFile out(outname + ".fastq", "");


        logfile().list(bed.size());
        logfile().list(bed.length());
        logfile().list(bed[0][0]);
        logfile().list(bed[0][1]);
        logfile().list(bed[1][0]);


        //out.writeln("the first three bases");
        //out.writeln(reference(0,0));
        //out.writeln(reference(0,1));
        //out.writeln(reference(0,2));
        //out.writeln(reference.view(0,0,2)); // this works


        for (const auto& b: bed) {
            b.refID();
            //b.pos(); // I don't think that exists, because b is a GenomeWindow instead of a GenomePosition

            b.size(); // we can check if the length of the bed thing is, so we can exclude all the ones that are not one

            // genotype position instead of genotype window?
            //genometools::TGenomePosition center(b.refID(), b.fromOnChr());
            //out.writeln("check if center works");
            //out.writeln(reference[center]); // works

            //out.writeln(reference(b.refID(), b.fromOnChr())); // this works and gives me the correct starting base
            //out.writeln(reference.view(b)); // this also works but with spaces

            if (flank > b.fromOnChr()) {
                continue;
            }


            size_t start = b.fromOnChr() - flank;
            size_t count = 2*flank + 1;


            out.writeln("@", bed.name(b.refID()), "_", start, "-", start+count); // first thing we actually need for the fastq, but we need to add more info


            out.writeln(reference.view(b.refID(), start, count));
            
            //out.writeln(reference(0,0));
            out.writeln("+");
            out.writeln(std::string(2*flank+1, 'F'));

        }

        out.close();

    }
}

