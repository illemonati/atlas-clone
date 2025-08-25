#ifndef TFASTATOFASTQ_H_
#define TFASTATOFASTQ_H_

#include "TAlignmentTraverser.h"
#include "genometools/TBed.h"


namespace PopulationTools {

class TBedToFastq {
public:
	TBedToFastq(){};
	void run();
};

class TBamToBed final {
private: 
	BAM::TAlignmentTraverser _alnTraverser;
	genometools::TBedWithInfo<std::string> _outBed;
	size_t _numPosNotAligned;
	
	void _traverserAlignments();

public:
	TBamToBed(){};
	void run();
};


class TPositionBasedLiftOver{
public:
	TPositionBasedLiftOver(){};
	void run();

};



} // namespace PopulationTools

#endif
