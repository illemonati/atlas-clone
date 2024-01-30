
#ifndef GENOMETASKS_TCONTEXTERRORS_H_
#define GENOMETASKS_TCONTEXTERRORS_H_

#include <string>

#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"

namespace GenomeTasks{

class TContextErrors final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
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
