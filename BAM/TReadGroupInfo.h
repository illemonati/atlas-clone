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
#include "stringFunctions.h"
#include "TParameters.h"
#include "TLog.h"
#include "TReadGroups.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"

namespace BAM {

namespace RGInfo{

//------------------------------------------------
// TInfoValue
//------------------------------------------------
//TODO: find a better way
enum class InfoType {RGName=0, seqType, numCycles, fragmentLengthDistr, baseQualityDistr, mappingQualityDistr, softClipDistr, recal, rho, COUNT};

std::string infoType2ArgumentString(InfoType Info) {
	switch (Info) {
		case InfoType::RGName: return "readGroupname";
		case InfoType::seqType: return "seqType";
		case InfoType::numCycles: return "numCycles";
		case InfoType::fragmentLengthDistr: return "fragmentLengthDistr";
		case InfoType::baseQualityDistr: return "baseQualityDistr";
		case InfoType::mappingQualityDistr: return "mappingQualityDistr";
		case InfoType::softClipDistr: return "softClipDistr";
		case InfoType::recal: return "recal";
		case InfoType::rho: return "rho";
		case InfoType::COUNT: DEVERROR("InfoType::COUNT not meant to be used.");
	}
	DEVERROR("InfoType missing!");
}

std::string infoType2Description(InfoType Info) {
	switch (Info) {
		case InfoType::RGName: return "read group name";
		case InfoType::seqType: return "sequencing type";
		case InfoType::numCycles: return "num sequencing cycles";
		case InfoType::fragmentLengthDistr: return "fragment length distribution";
		case InfoType::baseQualityDistr: return "quality distribution";
		case InfoType::mappingQualityDistr: return "mapping quality distribution";
		case InfoType::softClipDistr: return "softclip distribution";
		case InfoType::recal: return "base quality score recalibration model";
		case InfoType::rho: return "base quality score recalibration rho";
		case InfoType::COUNT: DEVERROR("InfoType::COUNT not meant to be used.");
	}
	DEVERROR("InfoType missing!");
}

std::string infoType2DefaultValue(InfoType Info) {
	switch(Info){
		case InfoType::RGName: return "SimReadGroup";

		case InfoType::seqType: return "single";
		case InfoType::numCycles: return "150";
		case InfoType::fragmentLengthDistr: return "fixed(300)";
		case InfoType::baseQualityDistr: return "normal(30,10)[0,93]";
		case InfoType::mappingQualityDistr: return "normal(60,10)[1,255]";
		case InfoType::softClipDistr: return "-";
		case InfoType::recal: return "-";
		case InfoType::rho: return "-";
		case InfoType::COUNT: DEVERROR("InfoType::COUNT not meant to be used.");
	}
	DEVERROR("InfoType missing!");
}

//------------------------------------------------
// TReadGroupInfoEntry
//------------------------------------------------

class TReadGroupInfoEntry{
private:
	std::array<std::string, static_cast<size_t>(InfoType::COUNT)> _rgInfo;

public:
	TReadGroupInfoEntry(const std::string & Name){
		set(InfoType::RGName, Name);
	}

	void set(InfoType Type, const std::string & Value){
		if(Type == InfoType::COUNT){
			DEVERROR("InfoType::COUNT not meant to be used.");
		}
		_rgInfo[static_cast<size_t>(Type)] = Value;
	}

	std::string get(InfoType Type) const {
		if(Type == InfoType::COUNT){
			DEVERROR("InfoType::COUNT not meant to be used.");
		}
		return _rgInfo[static_cast<size_t>(Type)];
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
	static inline const std::string _RGInfoArgument = "readGroupInfo";
	static inline const std::string _numRGArgument = "numReadGroups";

	std::vector<TReadGroupInfoEntry> _info;

	void _setAllReadGroups(InfoType Info, const std::string & Val);
	void _readInfoFromCommandLine(InfoType Info);
	void _readInfoFromRGInfoFile(InfoType Info, const std::vector< std::vector<std::string> > & fileData, size_t col);
	void _setDefault(InfoType Info);
	void _readInfoPerReadGroup(InfoType Info);
	void _readInfoPerReadGroup(InfoType Info, const std::vector<std::string> & header, const std::vector< std::vector<std::string> > & fileData);
	void _initializeFromRGInfoFile();
	void _initializeFromCommandLine();
	void _matchReadGroups(BAM::TReadGroups & ReadGroups);

public:
	TReadGroupInfo();

	// either: read info from file and match with TReadGroups (used for analyzes)
	void readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups);

	// or: read info and fill TReadGroups (used for simulations)
	void readInfoAndCreateReadGroups(BAM::TReadGroups & ReadGroups);

	// getters
	std::vector<TReadGroupInfoEntry>::const_iterator cbegin(){
		return _info.cbegin();
	}
	std::vector<TReadGroupInfoEntry>::const_iterator cend(){
		return _info.cend();
	}

	size_t size() const {
		return _info.size();
	}

	const TReadGroupInfoEntry& operator[](size_t index) const {
		return _info[index];
	}


};

} //end namespace RGInfo

} //end namespace BAM


#endif /* BAM_TREADGROUPINFO_H_ */
