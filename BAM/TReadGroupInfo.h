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
#include "TStrongArray.h"
#include "stringFunctions.h"
#include "TParameters.h"
#include "TLog.h"
#include "TReadGroups.h"
#include "TError.h"
#include "TTask.h"
#include "nlohmann/json.hpp"

namespace BAM {

namespace RGInfo{

//------------------------------------------------
// TInfoValue
//------------------------------------------------
enum class InfoType {min=0, RGName=0, RGFrequency, seqType, cycles, fragmentLength, baseQuality, mappingQuality, softClipping, recal1, recal2, pmd, max};

//------------------------------------------------
// argument string, description and default for each info type
//------------------------------------------------

struct TInfo {
	std::string argument;
	std::string description;
	std::string defaults;
	TInfo() = default;
	TInfo(std::string_view Argument, std::string_view Description, std::string_view Defaults)
		: argument(std::move(Argument)), description(std::move(Description)), defaults(std::move(Defaults)) {}
};

inline const coretools::TStrongArray<TInfo, InfoType> infos = []() {
	coretools::TStrongArray<TInfo, InfoType> i;
	i[InfoType::RGName] = {"readGroup", "read group name", "SimReadGroup"};
	i[InfoType::RGFrequency] = {"frequency", "read group frequency", "1.0"};
	i[InfoType::seqType] = {"seqType", "sequencing type", "single"};
	i[InfoType::cycles] = {"seqCycles", "number of sequencing cycles", "100"};
	i[InfoType::fragmentLength] = {"fragmentLength", "fragment length distribution", "fixed(300)"};
	i[InfoType::baseQuality] = {"baseQuality", "base quality distribution", "normal(30,10)[0,93]"};
	i[InfoType::mappingQuality] = {"mappingQuality", "mapping quality distribution", "normal(60,10)[1,255]"};
	i[InfoType::softClipping] = {"softClipping", "soft clipping distribution", "-"};
	i[InfoType::recal1] = {"recal1", "base quality score recalibration model for 1st mate", "-"};
	i[InfoType::recal2] = {"recal2", "base quality score recalibration model for 2nd mate", "-"};
	i[InfoType::pmd] = {"pmd", "Postmortem damage model", "-"};
	return i;
}();


using nlohmann::ordered_json;
//------------------------------------------------
// TReadGroupInfoEntry
//------------------------------------------------
class TReadGroupInfoEntry{
private:
	std::map<InfoType, ordered_json> _info;

public:
	TReadGroupInfoEntry(const std::string & RgName){
		_info.insert_or_assign(InfoType::RGName, RgName);
	}

	void set(const InfoType Info, const ordered_json & Value){
		if(Value != "" && Value != "-"){
			_info.insert_or_assign(Info, Value);
		}
	}

	bool has(const InfoType Info) const {
		return _info.find(Info) != _info.end();
	}

	const ordered_json& get(const InfoType Info) const;

	std::string getString(const InfoType Info) const;

	std::string name() const { return getString(InfoType::RGName); }

	const ordered_json& operator[](const InfoType Info) const {
		return get(Info);
	}

	void write(coretools::TOutputFile & Out, const InfoType Info) const;
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
	ordered_json _json;
	std::string _filename;
	coretools::TStrongArray<bool, InfoType> _parsed;

	void _setAllReadGroups(InfoType Info, const std::string & Val);
	void _setDefault(InfoType Info);
	void _setFromCommandLine(InfoType Info);
	void _setFromRGInfoFile(InfoType Info);
	bool _readGroupExists(const std::string & Name);
	void _readFile(const std::string & Filename);
	void _createReadGroupInfoEntries(const BAM::TReadGroups & ReadGroups);

public:
	static inline const std::string _RGInfoArgument = "RGInfo";
	static inline const std::string _numRGArgument = "numReadGroups";

	TReadGroupInfo() = default;
	~TReadGroupInfo() = default;

	// either: read info from file and match with TReadGroups (used for analyzes)
	void readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups, const std::string & Filename = "");

	// or: read info and fill TReadGroups (used for simulations)
	BAM::TReadGroups readInfoAndCreateReadGroups();
	BAM::TReadGroups readInfoAndCreateReadGroups(const std::string & RgInfoFileName);

	//functions to parse certain info
	void parse(const InfoType Info);

	template <typename... Ts>
	void parse(const InfoType Info, Ts... others){
		parse(Info);
		parse(others...);
	}

	// getters
	std::vector<TReadGroupInfoEntry>::const_iterator cbegin(){
		return _info.cbegin();
	}
	std::vector<TReadGroupInfoEntry>::const_iterator cend(){
		return _info.cend();
	}

	size_t size() const noexcept {
		return _info.size();
	}

	bool hasInfo(const InfoType Info) const {
		return _parsed[Info];
	};

	const TReadGroupInfoEntry& operator[](uint16_t RGIndex) const {
		return _info[RGIndex];
	}

	bool has(size_t RGIndex, InfoType Info) const noexcept {
		return _info[RGIndex].has(Info);
	};

	const ordered_json& get(size_t RGIndex, const InfoType Info) const noexcept {
		return _info[RGIndex][Info];
	}

	const ordered_json& get(size_t RGIndex, const InfoType Info, const ordered_json& defValue) const noexcept {
		return has(RGIndex, Info) ? get(RGIndex, Info) : defValue;
	}

	std::string getString(size_t RGIndex, const InfoType Info) const noexcept {
		return get(RGIndex, Info).dump();
	}

	std::string getString(size_t RGIndex, const InfoType Info, const std::string & defValue) const noexcept {
		return has(RGIndex, Info) ? getString(RGIndex, Info) : defValue;
	}

	template <typename Container>
	void fillContainerPerReadGroup(Container & Vec, const InfoType Info) const{
		Vec.resize(size());
		for(size_t i = 0; i < size(); ++i){
			coretools::str::fillFromString(get(i, Info), Vec[i]);
		}
	};

	bool hasFile() const { return _json.size() > 0; };
	bool fileHasInfo(const InfoType Info) const;

	std::vector<std::string> getUnusedAttributesInFile();
	void warnAboutUnusedColumnsInFile();

	// setters
	void set(const uint16_t RGIndex, const InfoType Info, const ordered_json & Value);

	//writing
	void write(const std::string & Filename);
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
