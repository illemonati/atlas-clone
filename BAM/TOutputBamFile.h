
#ifndef BAM_TOUTPUTBAMFILE_H_
#define BAM_TOUTPUTBAMFILE_H_

#include "TAlignment.h"
#include "TBamFile.h"
#include "TReadGroups.h"
#include "api/BamWriter.h"
#include "coretools/Types/probability.h"
#include <string>

namespace BAM {

//------------------------------------------------
// TQualityAdjusterForWriting
// Manages the printing of quality scores when writing BAM files
//------------------------------------------------
class TQualityAdjusterForWriting{
private:
	static constexpr char _minBaseQuality = coretools::toChar(coretools::PhredInt::highest());
	static constexpr char _maxBaseQuality = coretools::toChar(coretools::PhredInt::highest());
	bool _adjust      = false;
	bool _binIllumina = false;
	char _minQual = _minBaseQuality;
	char _maxQual = _maxBaseQuality;

	char _adjustOneQuality(char Qual) const;

public:
	TQualityAdjusterForWriting();
	void binQualitiesIllumina();
	void limitRange(const coretools::TNumericRange<char> & Range);
	std::string rangeString();
	void adjustQualities(std::string & qualities) const;
};

class TOutputBamFile {
private:
	std::string _outputFilename;
	BamTools::BamWriter _bamWriter;
	const TReadGroups *_readGroups = nullptr;

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
