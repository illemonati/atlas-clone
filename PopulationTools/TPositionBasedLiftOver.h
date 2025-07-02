#ifndef TFASTATOFASTQ_H_
#define TFASTATOFASTQ_H_

#include "TBamTraverser.h"
#include "genometools/TBed.h"


namespace PopulationTools {

class TBedToFastq {
public:
	TBedToFastq(){};
	void run();
};

class TBamToBed : public GenomeTasks::TBamReadTraverser<GenomeTasks::ReadType::Parsed> {
private: 
	genometools::TBedWithInfo<std::string> _outBed;
	size_t _numPosNotAligned;
	
	void _handleAlignment(BAM::TAlignment& alignment) override;

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
