#include "TSimulatorBamFiles.h"

#include "TSamHeader.h"
#include "coretools/Main/TLog.h"

#include "genometools/GenomePositions/TChromosomes.h"

namespace Simulations {

using coretools::instances::logfile;
using coretools::str::toString;

TSimulatorBamFiles::TSimulatorBamFiles(size_t NumFiles, const std::string & Outname, std::vector<TReadSimulators> & ReadSimulators,
				       const genometools::TChromosomes &Chromosomes) {
	coretools::instances::logfile().startIndent("Preparing BAM files for output:");

	if (NumFiles < 1) DEVERROR("Can not open less than one BAM file!");
	if(ReadSimulators.size() > 1 && ReadSimulators.size() != NumFiles){
		DEVERROR("Number of read simulators does not match number of files!");
	}
	if (Chromosomes.size() < 1) UERROR("Can not open a BAM file without specified chromosomes!");

	//create header
	const BAM::TSamHeader header("1.6", "coordinate", "none");

	//read quality adjustment for wiriting
	logfile().startIndent("Output options:");
	BAM::TQualityAdjusterForWriting qualityAdjuster;
	logfile().endIndent();

	// open BAM files
	logfile().startIndent("Opening ", NumFiles, " BAM files:");
	_files.reserve(NumFiles);
	if (NumFiles == 1) {
		_createBamFile(Outname + ".bam", "Ind1", header, ReadSimulators.front().readGroups(), Chromosomes, qualityAdjuster);
	} else {
		//check if ReadSimulators are shared
		for (size_t i = 0; i < NumFiles; ++i) {
			std::string sampleName = toString(i + 1);
			std::string filename = Outname + "_ind" + sampleName + ".bam";
			if(ReadSimulators.size() == 1){
				_createBamFile(filename, sampleName, header, ReadSimulators.front().readGroups(), Chromosomes, qualityAdjuster);
			} else {
				_createBamFile(filename, sampleName, header, ReadSimulators[i].readGroups(), Chromosomes, qualityAdjuster);
			}
		}
		logfile().endIndent();
	}
	logfile().endIndent();
}

void TSimulatorBamFiles::_createBamFile(const std::string & Filename,
										const std::string & SampleName,
										const BAM::TSamHeader & Header,
										BAM::TReadGroups & ReadGroups,
										const genometools::TChromosomes & Chromosomes,
										const BAM::TQualityAdjusterForWriting & QualityAdjuster){
	logfile().list(Filename);
	// create read group object
	for (auto &rg : ReadGroups){
		rg.sample_SM = SampleName;
	}

	_files.emplace_back(Filename, Header, Chromosomes, ReadGroups, QualityAdjuster);

	//unset sample name
	//TODO: ugly hack, find way to copy TReadGroups
	for (auto &rg : ReadGroups) {
		rg.sample_SM = "";
	}
}

void TSimulatorBamFiles::close() {
	logfile().startIndent("Indexing BAM files:");
	_files.clear();
	logfile().endIndent();
}

BAM::TOutputBamFile &TSimulatorBamFiles::operator[](size_t i) {
	if (i >= _files.size()) UERROR("BAM file ", i, " does not exist!");
	return _files[i];
}

}
