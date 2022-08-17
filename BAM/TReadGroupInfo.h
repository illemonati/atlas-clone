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
#include "TPostMortemDamage.h"

namespace BAM {

namespace RGInfo{

//------------------------------------------------
// TInfoValue
//------------------------------------------------
enum class InfoType {min=0, RGName=0, seqType, numCycles, fragmentLengthDistr, baseQualityDistr, mappingQualityDistr, softClipDistr, recal, rho, max};

//------------------------------------------------
// TReadGroupInfoEntry
//------------------------------------------------

// Alternative:
using TReadGroupInfoEntry = coretools::TStrongArray<std::string, InfoType>;

/*
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
*/

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
	static inline const std::string _RGInfoArgument = "readGroupInfo";
	static inline const std::string _numRGArgument = "numReadGroups";

	TReadGroupInfo() = default;
	~TReadGroupInfo() = default;

	// either: read info from file and match with TReadGroups (used for analyzes)
	void readInfoAndMatchReadGroups(const BAM::TReadGroups & ReadGroups);

	// or: read info and fill TReadGroups (used for simulations)
	BAM::TReadGroups readInfoAndCreateReadGroups();

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

	const TReadGroupInfoEntry& operator[](size_t index) const {
		return _info[index];
	}
};

} //end namespace RGInfo

} //end namespace BAM


#endif /* BAM_TREADGROUPINFO_H_ */
