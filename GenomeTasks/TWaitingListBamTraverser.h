#ifndef TWAITINGLISTBAMTRAVERSER_H_
#define TWAITINGLISTBAMTRAVERSER_H_

#include "TAlignment.h"
#include "TGenome.h"
#include "TOutputBamFile.h"
#include "TParser.h"
#include "genometools/BED/TBed.h"

namespace GenomeTasks {

enum class AlignmentStatus {orphan, filterOut, ready};

struct TWaitingAlignment{
	BAM::TAlignment alignment;
	AlignmentStatus status;
	TWaitingAlignment(BAM::TAlignment Alignment, AlignmentStatus Status) : alignment(Alignment), status(Status) {}
};
inline bool operator<(const TWaitingAlignment& lhs, const TWaitingAlignment& rhs) {
	return lhs.alignment < rhs.alignment;
}

class TWaitingListBamTraverser {
public:
	using container = std::vector<TWaitingAlignment>;
	using iterator  = typename container::iterator;

private:
	bool _doMasking       = false;
	bool _considerRegions = false;
	genometools::TBed _mask;
	void _setMasks(const genometools::TChromosomes& Chromosomes);

protected:
	TGenome _genome;
	TParser _parser;

	BAM::TAlignmentList _blacklist; //used to keep track of filtered out mates
	container _waitingList;

	BAM::TOutputBamFile _outBam;

	size_t _maxDistanceBetweenMates;
	bool _recalibrate;
	bool _incorporatePMD;
	bool _keepOrphans;
	bool _removeSoftClippedBases;
	size_t _maxNumberOfSoftClippedBases;

	void _writeOrFilter(TWaitingAlignment& WAlignment);
	void _writeOrphan(TWaitingAlignment& WAlignment);
	void _writeAll();
	void _writeUpTo(const genometools::TGenomePosition & position);
	BAM::TAlignment _parseIntoNewAlignment();

	//pure virtual functions
	virtual void _handleMates(TWaitingAlignment &Mate) = 0;
	virtual void _handleSingle()                       = 0;
	virtual bool _alignmentCanBeWrittenUnchanged()     = 0;

public:
	TWaitingListBamTraverser(std::string_view OutName);
	void traverseBAM();
};

}

#endif
