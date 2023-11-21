#ifndef TBAMTRAVERSER_H_
#define TBAMTRAVERSER_H_

#include "TGenome.h"
#include "TParser.h"

namespace GenomeTasks {

class TFilteredBamTraverser  {
protected:
	TGenome _genome;
	void _traverseBAMPassedQC();
	virtual void _handleAlignment() = 0;

public:
	TFilteredBamTraverser();
};

class TParsedBamTraverser  {
protected:
	TGenome _genome;
	TParser _parser;

	void _traverseBAMPassedQC();
	virtual void _handleAlignment(BAM::TAlignment& alignment) = 0;

public:
	TParsedBamTraverser();
};

template<bool isParsed>
using TBamTraverser = std::conditional_t<isParsed, TParsedBamTraverser, TFilteredBamTraverser>;
}
#endif
