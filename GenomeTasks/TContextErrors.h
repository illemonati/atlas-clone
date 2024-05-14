#ifndef GENOMETASKS_TCONTEXTERRORS_H_
#define GENOMETASKS_TCONTEXTERRORS_H_

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"

namespace GenomeTasks{

class TContextErrors final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	struct TCount {
		size_t noError  = 0;
		size_t oneError = 0;
	};

	coretools::TStrongArray<
		coretools::TStrongArray<TCount, genometools::Base, coretools::index(genometools::Base::N) + 1>,
		genometools::Base>
	    _counts{};
	coretools::TStrongArray<size_t, genometools::Base> _major{};
	coretools::TStrongArray<size_t, genometools::Base, coretools::index(genometools::Base::N) + 1> _minor{};
	coretools::TOutputFile _out;
	size_t _depth = 10;
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}
public:
	TContextErrors();
	void run();
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
