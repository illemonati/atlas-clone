
#ifndef TSIMULATORBAMFILES_H_
#define TSIMULATORBAMFILES_H_

#include "TOutputBamFile.h"
#include "TReadSimulators.h"
#include <vector>

namespace Simulations {

class TSimulatorBamFiles {
private:
	std::vector<BAM::TOutputBamFile> _files;

	void _createBamFile(const std::string & Filename,
											const std::string & SampleName,
											const BAM::TSamHeader & Header,
											BAM::TReadGroups & ReadGroups,
											const genometools::TChromosomes & Chromosomes,
											const BAM::TQualityAdjusterForWriting & QualityAdjuster);

public:
	TSimulatorBamFiles(size_t NumFiles, const std::string & Outname, std::vector<TReadSimulators> & ReadSimulators, const genometools::TChromosomes &Chromosomes);

	BAM::TOutputBamFile &operator[](size_t i);
	void close();
};
}

#endif
