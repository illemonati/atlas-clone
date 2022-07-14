/*
 * TSimulatorReadGroupInfo.h
 *
 *  Created on: Jul 14, 2022
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TSIMULATORREADGROUPINFO_H_
#define SIMULATIONS_TSIMULATORREADGROUPINFO_H_

// TODO: turn into read group info also used by TGenome

#include <vector>
#include "stringFunctions.h"
#include "TParameters.h"
#include "TLog.h"
#include "TReadGroups.h"
#include "TPostMortemDamage.h"
#include "TSequencingErrorModels.h"

namespace Simulations {

namespace RGInfo{

//------------------------------------------------
// TInfoValue
//------------------------------------------------
//TODO: find a better way
enum class InfoType {RGName, seqType, numCycles, fragmentLengthDistr, baseQualityDistr, mappingQualityDistr, softClipDistr, COUNT};

std::string infoType2ArgumentString(InfoType Info) {
	switch (Info) {
		case InfoType::RGName: return "readGroupname";
		case InfoType::seqType: return "seqType";
		case InfoType::numCycles: return "numCycles";
		case InfoType::fragmentLengthDistr: return "fragmentLengthDistr";
		case InfoType::baseQualityDistr: return "baseQualityDistr";
		case InfoType::mappingQualityDistr: return "mappingQualityDistr";
		case InfoType::softClipDistr: return "softClipDistr";
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
		case InfoType::COUNT: DEVERROR("IfnoType::COUNT not meant to be used.");
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
		case InfoType::COUNT: DEVERROR("IfnoType::COUNT not meant to be used.");
	}
	DEVERROR("InfoType missing!");
}

//------------------------------------------------
// TSimulatorReadGroupInfoEntry
//------------------------------------------------

class TSimulatorReadGroupInfoEntry{
private:
	std::array<std::string, static_cast<size_t>(InfoType::COUNT)> _rgInfo;
public:

	void set(InfoType Type, const std::string & Value){
		if(Type == InfoType::COUNT){
			DEVERROR("IfnoType::COUNT not meant to be used.");
		}
		_rgInfo[static_cast<size_t>(Type)] = Value;
	}

	std::string get(InfoType Type) const {
		if(Type == InfoType::COUNT){
			DEVERROR("IfnoType::COUNT not meant to be used.");
		}
		return _rgInfo[static_cast<size_t>(Type)];
	}
};

//------------------------------------------------
// TSimulatorReadGroupInfo
//------------------------------------------------
class TSimulatorReadGroupInfo{
private:
	static inline const std::string _RGInfoArgument = "readGroupInfo";
	static inline const std::string _numRGArgument = "numReadGroups";
	BAM::TReadGroups _readGroups;
	std::vector<TSimulatorReadGroupInfoEntry> _info;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	GenotypeLikelihoods::SequencingError::TModels _recal;

	void _setAllReadGroups(InfoType Info, const std::string & Val);
	void _readInfoFromCommandLine(InfoType Info);
	void _readInfoFromFile(InfoType Info, const std::vector< std::vector<std::string> > & fileData, size_t col);
	void _setDefault(InfoType Info);
	void _readInfoPerReadGroup(InfoType Info);
	void _readInfoPerReadGroup(InfoType Info, const std::vector<std::string> & header, const std::vector< std::vector<std::string> > & fileData);
	void _initializeFromRGInfoFile();
	void _initializeFromCommandLine();

public:
	TSimulatorReadGroupInfo();

};

} //end namespace RGInfo

} //end namespace Simulations


#endif /* SIMULATIONS_TSIMULATORREADGROUPINFO_H_ */
