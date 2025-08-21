#ifndef BAM_TREAD_H_
#define BAM_TREAD_H_

#include "TCigar.h"
#include "api/BamAlignment.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace BAM{
struct TRead {
	BamTools::BamAlignment bamAlignment;
	TCigar cigar;
	genometools::TGenomePosition position;
	size_t readGroupID;
	bool QCFiltersPassed;

	const std::string &name() const { return bamAlignment.Name; }
	uint16_t mappingQuality() const noexcept { return bamAlignment.MapQuality; }
	bool isPaired() const { return bamAlignment.IsPaired(); }
	bool isProperPair() const { return bamAlignment.IsProperPair(); }
	bool isReverseStrand() const { return bamAlignment.IsReverseStrand(); }
	bool isDuplicate() const { return bamAlignment.IsDuplicate(); }
	bool isMapped() const { return bamAlignment.IsMapped(); }
	bool isFailedQC() const { return bamAlignment.IsFailedQC(); }
	bool isSecondary() const { return !bamAlignment.IsPrimaryAlignment(); }
	bool isSupplementary() const { return bamAlignment.IsSupplementary(); }
	bool isFirstMate() const { return bamAlignment.IsFirstMate(); }
	bool isSecondMate() const { return bamAlignment.IsSecondMate(); }
	const std::string querySequence() const {return bamAlignment.QueryBases;}

	size_t fragmentLength() const;
	bool isLongerThanFragment() const;
};
}

#endif
