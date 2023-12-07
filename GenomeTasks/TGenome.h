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
	TGenome(std::string_view Name, const BAM::TBamFilters& Filters, size_t i = -1);

	~TGenome();
	TGenome(TGenome&&) = default;
	TGenome(const TGenome&) = delete;
	TGenome& operator=(TGenome&&) = default;
	TGenome& operator=(const TGenome&) = delete;

	const BAM::TBamFile& bamFile() const noexcept {return _bamFile;}
	BAM::TBamFile& bamFile() noexcept {return _bamFile;}

	const BAM::RGInfo::TReadGroupInfo& rgInfo() const noexcept {return _rgInfo;}
	BAM::RGInfo::TReadGroupInfo& rgInfo() noexcept {return _rgInfo;}

	const GenotypeLikelihoods::TErrorModels& errorModels() const noexcept {return _errorModels;};

	const std::string& outputName() const noexcept {return _outputName;}
};

}; // namespace GenomeTasks

#endif /* GENOME_H_ */
