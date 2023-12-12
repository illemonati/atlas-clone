#ifndef TWAITINGLISTBAMTRAVERSER_H_
#define TWAITINGLISTBAMTRAVERSER_H_

#include "TAlignment.h"
#include "TGenome.h"
#include "TOutputBamFile.h"
#include "TParser.h"

namespace GenomeTasks {

class TAlignmentMergerEntry{
private:
	BAM::TAlignment* _alignment;
	bool _ready;

public:

	TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting);
	TAlignmentMergerEntry(TAlignmentMergerEntry && other);
	~TAlignmentMergerEntry();
	TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other);

	//getters
	bool ready() const { return _ready; }
	const BAM::TAlignment& alignment() const { return *_alignment; }
	BAM::TAlignment& alignment() { return *_alignment; }
	const std::string& name() const;

	//setters
	void makeReady() { _ready = true; }
	void setAsNonProperPair() ;
	bool operator<(const TAlignmentMergerEntry & other) const;
};

template<typename StorageType>
auto findInStorage(StorageType & Storage, const std::string & name){
	for(auto it=Storage.begin(); it!=Storage.end(); ++it){
		if(it->name() == name){
			return it;
		}
	}
	return Storage.end();
}

class TWaitingListBamTraverser {
public:
	using container = std::vector<TAlignmentMergerEntry>;
	using iterator  = typename container::iterator;

protected:
	TGenome _genome;
	TParser _parser;

	BAM::TAlignmentList _blacklist; //used to keep track of filtered out mates
	container _alignmentStorage;

	BAM::TOutputBamFile _outBam;

	size_t _maxDistanceBetweenMates;
	bool _recalibrate;
	bool _incorporatePMD;
	bool _keepOrphans;
	bool _removeSoftClippedBases;
	size_t _maxNumberOfSoftClippedBases;

	void _writeOrFilter(TAlignmentMergerEntry& Entry); 
	void _writeAll();
	void _writeUpTo(const genometools::TGenomePosition & position);
	BAM::TAlignment* _parseIntoNewAlignment();

	//pure virtual functions
	virtual void _handleMates(BAM::TAlignment &alignment, iterator mate) = 0;
	virtual void _handleSingle(BAM::TAlignment &alignment)               = 0;
	virtual bool _alignmentCanBeWrittenUnchanged()                       = 0;

public:
	TWaitingListBamTraverser(std::string_view OutName);
	void traverseBAM();
};

}

#endif
