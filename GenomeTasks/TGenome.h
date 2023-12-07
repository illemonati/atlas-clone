/*
 * TGenome.h
 */

#ifndef TGENOME_H_
#define TGENOME_H_

#include <string>

#include "TBamFile.h"
#include "TReadGroupInfo.h"
#include "TErrorModels.h"

namespace GenomeTasks {

class TGenome {
private:
	BAM::TBamFile _bamFile;
	std::string _outputName;
	BAM::RGInfo::TReadGroupInfo _rgInfo;
	GenotypeLikelihoods::TErrorModels _errorModels;

public:
	TGenome(const BAM::TBamFilters& Filters = {false});
	~TGenome();

	const BAM::TBamFile& bamFile() const noexcept {return _bamFile;}
	BAM::TBamFile& bamFile() noexcept {return _bamFile;}

	const BAM::RGInfo::TReadGroupInfo& rgInfo() const noexcept {return _rgInfo;}
	BAM::RGInfo::TReadGroupInfo& rgInfo() noexcept {return _rgInfo;}

	const GenotypeLikelihoods::TErrorModels& errorModels() const noexcept {return _errorModels;};

	const std::string& outputName() const noexcept {return _outputName;}
};

}; // namespace GenomeTasks

#endif /* GENOME_H_ */
