#ifndef TBAMTRAVERSER_H_
#define TBAMTRAVERSER_H_

#include "TGenome.h"
#include "TParser.h"

namespace GenomeTasks {

class TFilteredBamTraverser  {
protected:
	TGenome _genome{BAM::TBamFilters{true}};
	void _traverseBAMPassedQC();
	virtual void _handleAlignment() = 0;

public:
	virtual ~TFilteredBamTraverser() = default;
};

class TParsedBamTraverser  {
protected:
	TGenome _genome{BAM::TBamFilters{true}};
	TParser _parser;

	void _traverseBAMPassedQC();
	virtual void _handleAlignment(BAM::TAlignment& alignment) = 0;

public:
	virtual ~TParsedBamTraverser() = default;
};

enum class ReadType : bool {Filtered, Parsed};

template<ReadType Type>
using TBamReadTraverser = std::conditional_t<Type == ReadType::Parsed, TParsedBamTraverser, TFilteredBamTraverser>;
} // namespace GenomeTasks
#endif
