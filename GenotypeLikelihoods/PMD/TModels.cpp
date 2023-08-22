#include "TModels.h"

namespace GenotypeLikelihoods::PMD {
void TModels::initialize(size_t NReadGroups, std::string_view PMDString) {
	if (PMDString.empty() || PMDString == "-" || PMDString == "defaul") {
		for (size_t i = 0; i < NReadGroups; ++i) { _pModels.push_back(&_noPMD); }
	}
}


void TModels::initialize(BAM::RGInfo::TReadGroupInfo & RGInfo) {
		
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
		
}

void TModels::pool(const BAM::TReadGroupMap& rgMap) {
	if (!hasPMD()) UERROR("No point pooling models that do not recalibrate!");
	for (size_t rg = 0; rg < _pModels.size(); ++rg) {
		const auto pIndex = rgMap.pooledIndex(rg);
		if (pIndex != rg) {
			_pModels[rg] = _pModels[pIndex];
		}
	}
		
}
} // namespace GenotypeLikelihoods::PMD
