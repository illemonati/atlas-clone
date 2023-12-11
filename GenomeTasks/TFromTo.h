
#ifndef GENOMETASKS_TFROMTO_H_
#define GENOMETASKS_TFROMTO_H_

#include <string>

#include "coretools/Files/TOutputFile.h"

#include "TBamWindowTraverser.h"

namespace GenomeTasks{

//----------------------------------------
// TDepthWriter
//----------------------------------------
class TFromTo final : public TBamWindowTraverser<WindowType::SingleBam> {
private:
	coretools::TOutputFile _out;
	void _handleWindow(GenotypeLikelihoods::TWindow& window) override;
	void _startChromosome(const genometools::TChromosome&) override {}
	void _endChromosome(const genometools::TChromosome&) override {}
public:
	TFromTo();
	void run();
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
