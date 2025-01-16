#ifndef GENOMETASKS_TOVERLAPQUANTIFIER_H_
#define GENOMETASKS_TOVERLAPQUANTIFIER_H_

#include "TGenome.h"
#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {

class TOverlapQuantifier {
private:
	TGenome _genome;
	std::vector<TWaitingAlignment> _alignmentStorage;

public:
	TOverlapQuantifier();
	void run();
};

}

#endif
