#ifndef GENOMETASKS_TTRANSITIONTABLER_H_
#define GENOMETASKS_TTRANSITIONTABLER_H_

#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {

class TTransitionTabler final
	: public TWaitingListBamTraverser {
private:
	using Rho         = coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base>, genometools::Base>;
	using Transitions = coretools::TStrongArray<std::vector<Rho>, BAM::Strand>;

	Transitions _transitions;

	coretools::TOutputFile _outTransitions;
	coretools::TOutputFile _outTransitionsRel;
	coretools::TOutputFile _outTransitionsPsi;
	coretools::TOutputFile _outTransitionsRho;

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;

	void _writeTransitions();

public:
	TTransitionTabler();
	void run();
};

}
#endif
