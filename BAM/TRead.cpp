#include "TRead.h"

namespace BAM {

size_t TRead::fragmentLength() const{
	if(bamAlignment.IsProperPair()){
		const size_t inserted = abs(bamAlignment.InsertSize) + cigar.lengthInserted();
		const size_t deleted  = cigar.lengthDeleted();
		if (inserted < deleted) return 0;
		return inserted - deleted;
	} else {
		return cigar.lengthSequenced();
	}
}

bool TRead::isLongerThanFragment() const {
	return bamAlignment.IsProperPair() && bamAlignment.InsertSize < static_cast<int>(cigar.lengthAligned());
}

} // namespace BAM
