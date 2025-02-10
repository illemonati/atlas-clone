#ifndef GENOMETASKS_TTRANSITIONTABLER_H_
#define GENOMETASKS_TTRANSITIONTABLER_H_

#include "TWaitingListBamTraverser.h"

namespace GenomeTasks {

class TTransitionTabler final
	: public TWaitingListBamTraverser {
private:
	using Rho         = coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base, 5>, genometools::Base, 5>;
	using Transitions = coretools::TStrongArray<std::vector<Rho>, BAM::Strand>;

	std::array<Transitions, 3> _transitions{};

	std::array<coretools::TOutputFile, 3> _outTransitions;
	std::array<coretools::TOutputFile, 3> _outTransitionsRel;
	std::array<coretools::TOutputFile, 3> _outTransitionsRho;

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;

	void _writeTransitions();

public:
	TTransitionTabler();
	void run();
};

}
#endif
