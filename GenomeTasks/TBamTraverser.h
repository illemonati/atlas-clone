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
	TFilteredBamTraverser() : _genome(BAM::TBamFilters{true}) {};
};

class TParsedBamTraverser  {
protected:
	TGenome _genome;
	TParser _parser;

	void _traverseBAMPassedQC();
	virtual void _handleAlignment(BAM::TAlignment& alignment) = 0;

public:
	TParsedBamTraverser() : _genome(BAM::TBamFilters{true}) {};
};

enum class ReadType : bool {Filtered, Parsed};

template<ReadType Type>
using TBamReadTraverser = std::conditional_t<Type == ReadType::Parsed, TParsedBamTraverser, TFilteredBamTraverser>;
} // namespace GenomeTasks
#endif
