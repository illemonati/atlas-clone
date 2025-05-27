#include "TModels.h"
#include "TReadGroupInfo.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"

namespace GenotypeLikelihoods::PMD {
using coretools::instances::logfile;

void TModels::initialize(size_t NReadGroups, const BAM::TReadGroupMap &rgMap) {
	_withPMD.reserve(NReadGroups);
	for (size_t i = 0; i < NReadGroups; ++i) {
		_withPMD.emplace_back();
		_pModels.push_back(&_withPMD.back());
	}
	_pool(rgMap);
}

void TModels::initialize(BAM::RGInfo::TReadGroupInfo & RgInfo) {
	using BAM::RGInfo::InfoType;
	std::vector<int> iis(RgInfo.size(), -1);
	bool reFormat = false;
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		const auto &Info = RgInfo[rg];
		if (Info.has(InfoType::pmd)) {
			const auto &json = Info[InfoType::pmd];
			if (json.empty()) continue;
			if (json.is_string()) {
				reFormat = true;
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
	if (reFormat) addToRGInfo(RgInfo);
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
	for(size_t r = 0; r < _pModels.size(); ++r){
		RgInfo.set(r, BAM::RGInfo::InfoType::pmd, model(r).info());
	}
}

void TModels::log(size_t rgID) const {
	if (model(rgID).hasPMD()) {
		_pModels[rgID]->psi()->log();
	} else {
		logfile().list("No PMD.");
	}
}

void TModels::_pool(const BAM::TReadGroupMap& rgMap) {
	coretools::user_assert(hasPMD(), "No point pooling models that do not recalibrate!");
	for (size_t rg = 0; rg < _pModels.size(); ++rg) {
		const auto pIndex = rgMap.pooledIndex(rg);
		if (pIndex != rg) {
			_pModels[rg] = _pModels[pIndex];
		}
	}
		
}
} // namespace GenotypeLikelihoods::PMD
