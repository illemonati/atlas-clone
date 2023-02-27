/*
 * TSimulatorReadGroupInfo.h
 *
 *  Created on: Jul 14, 2022
 *      Author: phaentu
 */

#ifndef BAM_TREADGROUPINFO_H_
#define BAM_TREADGROUPINFO_H_

// TODO: turn into read group info also used by TGenome

#include <vector>
#include "coretools/Containers/TBitSet.h"
#include "coretools/Strings/toString.h"
#include "coretools/devtools.h"
#include "nlohmann/json.hpp"

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TTask.h"
#include "coretools/Strings/stringFunctions.h"

#include "TReadGroups.h"
namespace BAM {

namespace RGInfo{

//------------------------------------------------
// Info and functions to extract data
//------------------------------------------------

using TInfo = nlohmann::ordered_json;

//------------------------------------------------
// TInfoValue
//------------------------------------------------
enum class InfoType : size_t {min=0, RGName=0, RGFrequency, seqType, cycles, fragmentLength, baseQuality, mappingQuality, softClipping, recal, pmd, max};

//------------------------------------------------
// argument string, description and default for each info type
//------------------------------------------------

struct TInfoArgument {
	std::string argument;
	std::string description;
	std::string defaults;
	TInfoArgument() = default;
	TInfoArgument(std::string_view Argument, std::string_view Description, std::string_view Defaults)
		: argument(Argument), description(Description), defaults(Defaults) {}
};

inline const coretools::TStrongArray<TInfoArgument, InfoType> infos = []() {
	coretools::TStrongArray<TInfoArgument, InfoType> i;
	i[InfoType::RGName]         = {"readGroup", "read group name", "SimReadGroup"};
	i[InfoType::RGFrequency]    = {"frequency", "read group frequency", "1.0"};
	i[InfoType::seqType]        = {"seqType", "sequencing type", "single"};
	i[InfoType::cycles]         = {"seqCycles", "number of sequencing cycles", "100"};
	i[InfoType::fragmentLength] = {"fragmentLength", "fragment length distribution", "gamma(10,0.2)[30,200]"};
	i[InfoType::baseQuality]    = {"baseQuality", "base quality distribution", "normal(30,10)[0,93]"};
	i[InfoType::mappingQuality] = {"mappingQuality", "mapping quality distribution", "normal(60,10)[1,255]"};
	i[InfoType::softClipping]   = {"softClipping", "soft clipping distribution", "exponential(0.5)[0,20]"};
	i[InfoType::recal]          = {"recal", "base quality score recalibration model", "-"};
	i[InfoType::pmd]            = {"pmd", "Postmortem damage model", "-"};
	return i;
}();

//------------------------------------------------
// Predefined tags
//------------------------------------------------
namespace seqType{
	static constexpr std::string_view single = "single";
	static constexpr std::string_view paired = "paired";
}

//------------------------------------------------
// TReadGroupInfoEntry
//------------------------------------------------

class TReadGroupInfo;
void parse(TReadGroupInfo* rgi, InfoType i);

class TReadGroupInfoEntry{
private:
	TReadGroupInfo* _rgi;
	coretools::TStrongArray<TInfo, InfoType> _info{};

public:
	TReadGroupInfoEntry(TReadGroupInfo* rgi, std::string_view RgName) : _rgi(rgi) {
		_info[InfoType::RGName] = RgName;
	}

	bool has(InfoType Info) const {
		parse(_rgi, Info);
		return (!_info[Info].is_null() && _info[Info] != "-");
	}

	std::string getString(InfoType Info) const {
		parse(_rgi, Info);
		return coretools::str::toString(_info[Info]);
	};

	std::string name() const { return getString(InfoType::RGName); }

	const TInfo& operator[](InfoType Info) const {
		parse(_rgi, Info);
		return _info[Info];
	}
	
	TInfo& operator[](InfoType Info) noexcept {
		return _info[Info];
	}
};

//------------------------------------------------
// TReadGroupInfo
// Can be initialized from file, command line or default values in this order:
// 1) value provided on command line (could be a file name)
// 2) value provided in RG info file (error if RG is missing in file!)
// 3) default value
//
//------------------------------------------------
class TReadGroupInfo{
private:
	std::vector<TReadGroupInfoEntry> _info;
	TInfo _json;
	std::string _filename;
	coretools::TStrongBitSet<InfoType> _parsed;

	void _setAllReadGroups(InfoType Info, std::string_view Val);
	void _setDefault(InfoType Info);
	void _setFromCommandLine(InfoType Info);
	void _setFromRGInfoFile(InfoType Info);
	bool _readGroupExists(std::string_view Name);
	void _readFile(std::string_view Filename);
	void _createReadGroupInfoEntries(const BAM::TReadGroups & ReadGroups);
	void _parse(InfoType Info);

public:
	static constexpr std::string_view RGInfoArgument = "RGInfo";
	static constexpr std::string_view numRGArgument = "numReadGroups";

	TReadGroupInfo() = default;
	TReadGroupInfo(const BAM::TReadGroups & ReadGroups);
	TReadGroupInfo(const BAM::TReadGroups & ReadGroups, std::string_view Filename);

	// or: read info and fill TReadGroups (used for simulations)
	BAM::TReadGroups createReadGroups(std::string_view RgInfoFileName = "");

	// getters
	auto begin() const noexcept {
		return _info.begin();
	}
	auto end() const noexcept {
		return _info.end();
	}

	size_t size() const noexcept {
		return _info.size();
	}

	const TReadGroupInfoEntry& operator[](uint16_t RGIndex) const {
		assert(RGIndex < _info.size());
		return _info[RGIndex];
	}

	TReadGroupInfoEntry& operator[](uint16_t RGIndex) {
		assert(RGIndex < _info.size());
		return _info[RGIndex];
	}

	bool has(size_t RGIndex, InfoType Info) const noexcept {
		assert(RGIndex < _info.size());
		return _info[RGIndex].has(Info);
	};

	const TInfo& get(size_t RGIndex, const InfoType Info) const noexcept {
		assert(RGIndex < _info.size());
		return _info[RGIndex][Info];
	}

	const TInfo& get(size_t RGIndex, const InfoType Info, const TInfo& defValue) const noexcept {
		return has(RGIndex, Info) ? get(RGIndex, Info) : defValue;
	}

	std::string getString(size_t RGIndex, const InfoType Info) const noexcept {
		return _info[RGIndex].getString(Info);
	}

	std::string getString(size_t RGIndex, const InfoType Info, std::string_view defValue) const noexcept {
		return has(RGIndex, Info) ? getString(RGIndex, Info) : std::string{defValue};
	}

	template <typename Container>
	void fillContainerPerReadGroup(Container & Vec, const InfoType Info) const{
		Vec.resize(size());
		for(size_t i = 0; i < size(); ++i){
			coretools::str::fromString(get(i, Info).get<std::string_view>(), Vec[i]);
		}
	};

	bool hasFile() const { return _json.size() > 0; };
	bool fileHasInfo(const InfoType Info) const;

	std::vector<std::string> getUnusedAttributesInFile();
	void warnAboutUnusedColumnsInFile();

	// setters
	void set(const uint16_t RGIndex, const InfoType Info, const TInfo & Value);

	//writing
	void write(std::string_view Filename);

	// preparse arguments
	void parse(InfoType Info){ _parse(Info); }
	template <typename... InfoTypes>
	void parse(InfoType Info, InfoTypes... FurtherInfos){
		_parse(Info);
		parse(FurtherInfos...);
	}

	friend void parse(TReadGroupInfo *ReadGroupInfo, InfoType Info) { ReadGroupInfo->_parse(Info); }
};


//-------------------------------------------
// TTask_testReadGroupInfo
//-------------------------------------------
class TTask_testReadGroupInfo : public coretools::TTask{
public:
	TTask_testReadGroupInfo(){ _explanation = "Testing JSON stuff"; };

	void run();
};


} //end namespace RGInfo

} //end namespace BAM


#endif /* BAM_TREADGROUPINFO_H_ */
