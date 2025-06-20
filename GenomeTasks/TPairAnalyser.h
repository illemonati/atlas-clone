#ifndef GENOMETASKS_TPAIRANALYSER_H_
#define GENOMETASKS_TPAIRANALYSER_H_

#include "TWaitingListBamTraverser.h"
#include "coretools/Files/TOutputFile.h"

namespace GenomeTasks {

class TPairAnalyser final
	: public TWaitingListBamTraverser {
private:
	coretools::TOutputFile _out;

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &) override {}
	void _handleOrphan(TWaitingAlignment &) override {}

public:
	TPairAnalyser();
	void run() { traverseBAM(); }
};
}
#endif
