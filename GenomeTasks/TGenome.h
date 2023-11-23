/*
 * TGenome.h
 */

#ifndef TGENOME_H_
#define TGENOME_H_

#include <string>

#include "TBamFile.h"
#include "TReadGroupInfo.h"

namespace GenomeTasks {

//---------------------------------------------------------------
// TGenome_basic
// A base class without filters and genotype likelihoods
//---------------------------------------------------------------
class TGenome {
private:
	BAM::TBamFile _bamFile;
	std::string _outputName;
	BAM::RGInfo::TReadGroupInfo _rgInfo;

public:
	TGenome();
	~TGenome();

	const BAM::TBamFile& bamFile() const noexcept {return _bamFile;}
	BAM::TBamFile& bamFile() noexcept {return _bamFile;}

	const BAM::RGInfo::TReadGroupInfo& rgInfo() const noexcept {return _rgInfo;}
	BAM::RGInfo::TReadGroupInfo& rgInfo() noexcept {return _rgInfo;}

	const std::string& outputName() const noexcept {return _outputName;}
};

}; // namespace GenomeTasks

#endif /* GENOME_H_ */
