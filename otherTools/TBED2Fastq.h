
#ifndef BED2FASTQ_H
#define BED2FASTQ_H

#include "TParameters.h"
#include "TLog.h"
#include "TFastaReader.h"
#include "TBed.h"
#include "TLineWriter.h"

namespace Other{

class TBED2Fastq {
private:
	genometools::TBed _bed;
	genometools::TFastaReader _fasta;
	coretools::TLineWriter _out;

	
	
public:
	TBED2Fastq(std::string BedFileName, std::string FastaFileName, std::string OutFileName, int FlankingLength);

	void run();
};



struct TBED2FastqWriter {
	void run() {
		using coretools::instances::parameters;
		using coretools::instances::logfile;
		
        std::string bed = parameters().get<std::string>("bed");
		std::string fasta = parameters().get<std::string>("fasta");
		int flankingLength = parameters().get<int>("flank");
		std::string outName;
		if(parameters().exists("out")){
			outName = parameters().get<std::string>("out");
		} else {
			outName = coretools::str::readBeforeLast(bed, ".");
			outName += ".fastq";
		}

		// create object
		TBED2Fastq(bed, fasta, outName, flankingLength);		

	}
};

} // end namespcae Other
#endif // BED2FASTQ_H