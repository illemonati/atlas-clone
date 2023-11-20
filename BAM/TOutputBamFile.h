
#ifndef BAM_TOUTPUTBAMFILE_H_
#define BAM_TOUTPUTBAMFILE_H_

#include "TAlignment.h"
#include "TBamFile.h"
#include "TReadGroups.h"
#include "api/BamWriter.h"
#include <string>
namespace BAM {

//------------------------------------------------
// TQualityAdjusterForWriting
// Manages the printing of quality scores when writing BAM files
//------------------------------------------------
class TQualityAdjusterForWriting{
private:
	bool _adjust      = false;
	bool _binIllumina = false;
	genometools::BaseQuality _minQual{genometools::BaseQuality::min()};
	genometools::BaseQuality _maxQual {genometools::BaseQuality::max()};

	char _adjustOneQuality(genometools::BaseQuality qual) const;

public:
	TQualityAdjusterForWriting();
	void binQualitiesIllumina();
	void limitRange(const genometools::BaseQuality & min, const genometools::BaseQuality & max);
	void limitRange(const coretools::TNumericRange<uint8_t> & Range);
	std::string rangeString();
	void adjustQualities(std::string & qualities) const;
};

class TOutputBamFile {
private:
	std::string _outputFilename;
	BamTools::BamWriter _bamWriter;
	const TReadGroups *_readGroups = nullptr;

	// std::multiset<TAlignment, std::less<>> _futureAlignments;
	std::vector<TAlignment> _futureAlignments;

	// quality output transformations
	TQualityAdjusterForWriting _qualityAdjuster;

	void _writeAlignment(const TAlignment &alignment);

public:
	TOutputBamFile(std::string Filename, const TBamFile &Original);
	TOutputBamFile(std::string Filename, const TSamHeader &Header, const genometools::TChromosomes &Chromosomes,
				   const TReadGroups &ReadGroups);
	TOutputBamFile(std::string Filename, const TSamHeader &Header, const genometools::TChromosomes &Chromosomes,
				   const TReadGroups &ReadGroups, const TQualityAdjusterForWriting &QualityAdjuster);
	~TOutputBamFile();

	TOutputBamFile(const TOutputBamFile &)                = delete;
	TOutputBamFile &operator=(const TOutputBamFile &)     = delete;
	TOutputBamFile(TOutputBamFile &&) noexcept            = default;
	TOutputBamFile &operator=(TOutputBamFile &&) noexcept = default;

	void writeAlignment(const TAlignment &alignment);
	void writeAlignment(BamTools::BamAlignment &alignment);
	void writeAlignmentLater(const TAlignment &alignment);
};

} // namespace BAM
#endif
