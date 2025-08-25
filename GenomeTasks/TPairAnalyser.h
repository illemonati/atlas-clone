#ifndef GENOMETASKS_TPAIRANALYSER_H_
#define GENOMETASKS_TPAIRANALYSER_H_

#include "TWaitingListBamTraverser.h"
#include "coretools/Files/TOutputFile.h"

namespace GenomeTasks {

class TPairAnalyser final
	: public BAM::TWaitingListBamTraverser {
private:
	coretools::TOutputFile _out;

	void _handleMates(BAM::TWaitingAlignment &lhs, BAM::TWaitingAlignment &rhs) override;
	void _handleSingle(BAM::TWaitingAlignment &) override {}
	void _handleOrphan(BAM::TWaitingAlignment &) override {}

public:
	TPairAnalyser();
	void run() { traverseBAM(); }
};
}
#endif
