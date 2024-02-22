/*
 * TModels.cpp
 *
 *  Created on: May 14, 2020
 *      Author: phaentu
 */

#include "SequencingError/TModels.h"

namespace GenotypeLikelihoods {
namespace SequencingError {

using BAM::Mate;

//--------------------------------------------------------------------
// TModels
//--------------------------------------------------------------------

namespace impl {

std::pair<std::string_view, std::string_view> epsRho(std::string_view s) {
	// Format: intercept[];cov1:function1[];cov2:function2[];...;rho[[]]
	const auto rBegin = s.find("rho");
	if (rBegin == std::string::npos) {
		// no rho definition
		return std::make_pair(s, "default");
	}
	return std::make_pair(s.substr(0, rBegin-1), s.substr(rBegin + 3, s.size()));
}

} // namespace impl

void TModels::_pool(const BAM::TReadGroupMap& rgMap) {
	if (!recalibrates()) UERROR("No point pooling models that do not recalibrate!");
	for (size_t rg = 0; rg < _pModels.size(); ++rg) {
		const auto pIndex = rgMap.pooledIndex(rg);
		if (pIndex != rg) {
			_pModels[rg].front() = _pModels[pIndex].front();
			_pModels[rg].back() = _pModels[pIndex].back();
		}
	}
}

void TModels::initialize(size_t NReadGroups, std::string_view RecalString, const BAM::TReadGroupMap &rgMap) {
	_withRecal.reserve(NReadGroups * 2); // 2 mates per readgroup
	_pModels.reserve(NReadGroups);
	for (size_t i = 0; i < NReadGroups; ++i) {
		auto &first  = _withRecal.emplace_back(RecalString);
		auto &second = _withRecal.emplace_back(RecalString);
		_pModels.push_back(RGModels({&first, &second}));
	}
	_pool(rgMap);
}

void TModels::initialize(const BAM::RGInfo::TReadGroupInfo &RgInfo) {
	using BAM::RGInfo::InfoType;
	std::vector<coretools::TStrongArray<int, Mate>> iis(RgInfo.size(), {{-1, -1}});
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		const auto &Info = RgInfo[rg];

		// check if recal is provided
		if (Info.has(InfoType::recal)) {
			const auto &json = Info[InfoType::recal];
			if (json.contains("Mate1")) {
				const auto &info = json["Mate1"];
				if (info.empty()) {
					iis[rg].front() = -1;
				} else if (info.is_string()) {
					const auto sinfo = info.get<std::string_view>();
					if (sinfo.empty() || sinfo == "-" || sinfo == "default") {
						iis[rg].front() = -1;
					} else {
						auto [recal, rho] = impl::epsRho(info.get<std::string_view>());
						_withRecal.emplace_back(recal, rho);
						iis[rg].front() = _withRecal.size() - 1;
					}
				} else {
					_withRecal.emplace_back(info);
					iis[rg].front() = _withRecal.size() - 1;
				}
				if (json.contains("Mate2")) {
					const auto &info = json["Mate2"];
					if (info.empty()) {
						iis[rg].back() = -1;
					} else if (info.is_string()) {
						const auto sinfo = info.get<std::string_view>();
						if (sinfo.empty() || sinfo == "-" || sinfo == "default") {
							iis[rg].back() = -1;
						} else {
							auto [recal, rho] = impl::epsRho(info.get<std::string_view>());
							_withRecal.emplace_back(recal, rho);
							iis[rg].back() = _withRecal.size() - 1;
						}
					} else {
						_withRecal.emplace_back(info);
						iis[rg].back() = _withRecal.size() - 1;
					}
				} else {
					iis[rg].back() = -1;
				}
			} else {
				iis[rg].back() = -1; // no second mate
				if (json.empty()) {
					iis[rg].front() = -1;
				} else if (json.is_string()) {
					const auto sinfo = json.get<std::string_view>();
					if (sinfo.empty() || sinfo == "-" || sinfo == "default") {
						iis[rg].front() = -1;
					} else {
						auto [recal, rho] = impl::epsRho(json.get<std::string_view>());
						_withRecal.emplace_back(recal, rho);
						iis[rg].front() = _withRecal.size() - 1;
					}
				} else {
					_withRecal.emplace_back(json);
					iis[rg].front() = _withRecal.size() - 1;
				}
			}
		} else {
			iis[rg].front() = -1;
			iis[rg].back() = -1;
		}
	}
	_pModels.resize(RgInfo.size());
	for (size_t rg = 0; rg < RgInfo.size(); ++rg) {
		for (Mate m = Mate::min; m < Mate::max; ++m) {
			if (iis[rg][m] == -1) _pModels[rg][m] = &_noRecal;
			else _pModels[rg][m] = &_withRecal[iis[rg][m]];
		}
	}
}

void TModels::addToRGInfo(BAM::RGInfo::TReadGroupInfo & RgInfo) const {
	for(size_t r = 0; r < _pModels.size(); ++r){
		const auto& rgModel = RGModel(r);
		if (RGModel(r).back()->recalibrates()) {
			BAM::RGInfo::TInfo info;
			info["Mate1"] = rgModel.front()->info();
			info["Mate2"] = rgModel.back()->info();
			RgInfo.set(r, BAM::RGInfo::InfoType::recal, info);
		} else {
			RgInfo.set(r, BAM::RGInfo::InfoType::recal, rgModel.front()->info());
		}
	}
}

} // namespace SequencingError
}; // namespace GenotypeLikelihoods
