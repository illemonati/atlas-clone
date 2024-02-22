#include "TModels.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods::PMD {
void TModels::initialize(size_t NReadGroups, std::string_view PMDString, const BAM::TReadGroupMap &rgMap) {
	_withPMD.reserve(NReadGroups);
	for (size_t i = 0; i < NReadGroups; ++i) {
		_withPMD.emplace_back(PMDString);
		_pModels.push_back(&_withPMD.back());
	}
	_pool(rgMap);
}

void TModels::initialize(const BAM::RGInfo::TReadGroupInfo & RgInfo) {
	using BAM::RGInfo::InfoType;
	std::vector<int> iis(RgInfo.size(), -1);
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		const auto &Info = RgInfo[rg];
		if (Info.has(InfoType::pmd)) {
			const auto &json = Info[InfoType::pmd];
			if (json.empty()) continue;
			if (json.is_string()) {
				const auto sinfo = json.get<std::string_view>();
				if (sinfo.empty() || sinfo == "-" || sinfo == "default") continue;
				_withPMD.emplace_back(sinfo);
				iis[rg] = _withPMD.size() - 1;
			} else {
				_withPMD.emplace_back(json);
				iis[rg] = _withPMD.size() - 1;
			}
		}
	}
	_pModels.resize(RgInfo.size());
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		if (iis[rg] == -1) _pModels[rg] = &_noPMD;
		else _pModels[rg] = &_withPMD[iis[rg]];
	}
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
	for(size_t r = 0; r < _pModels.size(); ++r){
		RgInfo.set(r, BAM::RGInfo::InfoType::pmd, model(r).info());
	}
}

void TModels::_pool(const BAM::TReadGroupMap& rgMap) {
	if (!hasPMD()) UERROR("No point pooling models that do not recalibrate!");
	for (size_t rg = 0; rg < _pModels.size(); ++rg) {
		const auto pIndex = rgMap.pooledIndex(rg);
		if (pIndex != rg) {
			_pModels[rg] = _pModels[pIndex];
		}
	}
		
}
} // namespace GenotypeLikelihoods::PMD
