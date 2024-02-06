#include "TReadGroupMap.h"
#include "coretools/Files/TInputFile.h"
#include "coretools/Main/TError.h"
#include <filesystem>
#include "coretools/Main/TLog.h"

namespace BAM{

using coretools::instances::logfile;
//---------------------------------------------------------------
//TReadGroupMap
//---------------------------------------------------------------
TReadGroupMap::TReadGroupMap(const TReadGroups & ReadGroups, std::string_view Type) {
	if(Type.empty()){
		_noPooling(ReadGroups);
	} else if (std::filesystem::exists(Type)) {
		_fromFile(ReadGroups, Type);
	} else if (Type == "all"){
		_poolAll(ReadGroups);
	} else {
		UERROR("Cannot understand readgroup map argument: '", Type, "'!");
	}
};

void TReadGroupMap::_resize(const TReadGroups & ReadGroups){
	_readGroupMap.resize(ReadGroups.size(), ReadGroupMapNotInitializedIndex);
	_reverseReadGroupMap.resize(ReadGroups.size());
};

void TReadGroupMap::_markAsInUse(size_t index){
	_readGroupMap[index] = index;
	_reverseReadGroupMap[index].push_back(index);
	_readGroupsInUse.push_back(index);
};

void TReadGroupMap::_noPooling(const TReadGroups & ReadGroups){
	logfile().list("Not pooling any readgroups");
	_resize(ReadGroups)	;
	for(size_t r = 0; r < ReadGroups.size(); ++r){
		_markAsInUse(r);
	}
};

void TReadGroupMap::_poolAll(const TReadGroups & ReadGroups){
	const auto pool = 0;
	logfile().list("Pool all readgroups with ", ReadGroups.getName(pool), ".");
	_resize(ReadGroups)	;
	_markAsInUse(pool);
	for(size_t rg = 1; rg < ReadGroups.size(); ++rg){
		_readGroupMap[rg] = pool;
		_reverseReadGroupMap[pool].push_back(rg);
	}
}

void TReadGroupMap::_fromFile(const TReadGroups & ReadGroups, std::string_view filename){
	//set all values to no-initialized
	_resize(ReadGroups);

	//read read groups and their expected lengths
	logfile().listFlush("Reading read groups to be pooled from file '", filename, "' ...");
	coretools::TInputFile in(filename, coretools::FileType::Header, "\t", "//");

	std::vector< std::vector<std::string> > readGroupsToMerge;

	//parse file and fill vectors
	std::string readGroup;
	bool pooledAtLeastOneRG = false;

	for (; !in.empty(); in.popFront()) {
		//ignore if read group does not exist
		const auto rgName = in.get("readGroup");
		if(ReadGroups.readGroupExists(rgName)){
			//get read group index
			const auto poolWith = in.get("poolWith");
			size_t rg = ReadGroups.getId(rgName);
			size_t pool = ReadGroups.getId(poolWith);

			//check if rg to pool with is pooled itself
			if(_readGroupMap[pool] == ReadGroupMapNotInitializedIndex){
				_markAsInUse(pool);
			} else if(_readGroupMap[pool] != pool){
				UERROR("Read group '", poolWith, "' can not be pooled and pooled with at the same time!");
			}

			//check if rg can be pooled: allow self-pooling
			if(rg != pool){
				if(_readGroupMap[rg] == rg){
					UERROR("Read group '", rgName, "' can not be pooled and pool with at the same time!");
				} else if(_readGroupMap[rg] != ReadGroupMapNotInitializedIndex){
					UERROR("Read group '", rgName, "' is pooled multiple times in file '", filename, "'!");
				}

				//pool!
				_readGroupMap[rg] = pool;
				_reverseReadGroupMap[pool].push_back(rg);
				pooledAtLeastOneRG = true;
			}
		}
	}

	//mark all read groups not in file as pooled with itself
	for(size_t r = 0; r < _readGroupMap.size(); ++r){
		if(_readGroupMap[r] == ReadGroupMapNotInitializedIndex){
			_markAsInUse(r);
		}
	}
	logfile().done();

	//report
	if(pooledAtLeastOneRG){
		logfile().startIndent("The read groups will be pooled for parameter estimation as follows:");
		for(size_t r = 0; r < _readGroupMap.size(); ++r){
			if(_reverseReadGroupMap[r].size() > 0){
				logfile().startIndent(ReadGroups.getName(r) + ":");
				for(auto& s: _reverseReadGroupMap[r]){
					logfile().list(ReadGroups.getName(s));
				}
				logfile().endIndent();
			}
		}
	} else {
		logfile().warning("No read groups are pooled! Are you using the correct pooling file?");
	}
};

}
