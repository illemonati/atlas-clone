#ifndef TREADTRAVERSER_H_
#define TREADTRAVERSER_H_

#include "TBamFile.h"
#include "TErrorModels.h"
#include "TReadGroupInfo.h"

namespace BAM {
	
class TReadTraverser {
	BAM::TBamFile _bamFile;
	std::string _outputName;
	BAM::RGInfo::TReadGroupInfo _rgInfo;
	GenotypeLikelihoods::TErrorModels _errorModels;

	bool _eor = false;

public:
	TReadTraverser(bool EnableFilters=true);
	TReadTraverser(std::string_view Name, bool EnableFilters, size_t i);

	~TReadTraverser();
	TReadTraverser(TReadTraverser&&) = default;
	TReadTraverser(const TReadTraverser&) = delete;
	TReadTraverser& operator=(TReadTraverser&&) = default;
	TReadTraverser& operator=(const TReadTraverser&) = delete;

	const BAM::TBamFile &bamFile() const noexcept { return _bamFile; }
	BAM::TBamFile &bamFile() noexcept { return _bamFile; }

	const BAM::RGInfo::TReadGroupInfo &rgInfo() const noexcept { return _rgInfo; }
	BAM::RGInfo::TReadGroupInfo &rgInfo() noexcept { return _rgInfo; }

	const GenotypeLikelihoods::TErrorModels &errorModels() const noexcept { return _errorModels; };

	const std::string &outputName() const noexcept { return _outputName; }

	const genometools::TChromosomes &chromosomes() const noexcept { return _bamFile.chromosomes(); }
	const genometools::TChromosome &curChr() const noexcept { return _bamFile.curChromosome(); }

	const BAM::TRead &read() const noexcept { return bamFile().curRead(); }

	void nextRead();
	bool endOfReads();
};
} // namespace GenomeTasks

#endif
