#ifndef TWAITINGLISTBAMTRAVERSER_H_
#define TWAITINGLISTBAMTRAVERSER_H_

#include <memory>

#include "TAlignment.h"
#include "TOutputBamFile.h"
#include "TParser.h"
#include "TBamReadMask.h"
#include "TReadTraverser.h"

namespace BAM {

enum class AlignmentStatus {waiting, orphan, filterOut, ready};

struct TWaitingAlignment{
	BAM::TAlignment alignment;
	AlignmentStatus status = AlignmentStatus::waiting;
	TWaitingAlignment() = default;
	TWaitingAlignment(BAM::TAlignment Alignment, AlignmentStatus Status = AlignmentStatus::waiting)
		: alignment(std::move(Alignment)), status(Status) {}
};

class TWaitingListBamTraverser {
private:
	TBamReadMask _masks;

protected:
	TReadTraverser _genome{true};
	TParser _parser;

	BAM::TAlignmentList _blacklist; //used to keep track of filtered out mates
	std::vector<TWaitingAlignment> _waitingList;

	std::unique_ptr<BAM::TOutputBamFile> _outBam;

	size_t _maxDistanceBetweenMates;
	bool _recalibrate;
	bool _incorporatePMD;
	bool _keepOrphans;
	bool _removeSoftClippedBases;
	size_t _maxNumberOfSoftClippedBases;

	void _writeOrFilter(TWaitingAlignment& WAlignment);
	void _writeAll();
	void _writeUpTo(const genometools::TGenomePosition & position);
	TWaitingAlignment _nextAlignment();

	//pure virtual functions
	virtual void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) = 0;
	virtual void _handleSingle(TWaitingAlignment &lhs)                        = 0;
	virtual void _handleOrphan(TWaitingAlignment &lhs)                        = 0;
	bool _alignmentCanBeWrittenUnchanged();

public:
	TWaitingListBamTraverser(std::string_view OutName="");
	void traverseBAM();
};

}

#endif
