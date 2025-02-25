#ifndef GENOMETASKS_TTRANSITIONTABLER_H_
#define GENOMETASKS_TTRANSITIONTABLER_H_

#include "TAlignment.h"
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

	void _handleAlignments(const BAM::TAlignment& first, const BAM::TAlignment& second);
	void _handleFwd(std::vector<Rho>& transAll, std::vector<Rho>& transMate, const BAM::TAlignment& aln, size_t addL);
	void _handleRev(std::vector<Rho>& transAll, std::vector<Rho>& transMate, const BAM::TAlignment& aln, size_t addL);

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;

	void _writeTransitions();

public:
	TTransitionTabler();
	void run();
};

}
#endif
