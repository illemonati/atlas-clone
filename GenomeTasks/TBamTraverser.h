#ifndef TBAMTRAVERSER_H_
#define TBAMTRAVERSER_H_

#include "TParser.h"
#include "TReadTraverser.h"

namespace GenomeTasks {

class TParsedBamTraverser  {
protected:
	TReadTraverser _genome{true};
	TParser _parser;

	void _traverseBAMPassedQC();
	virtual void _handleAlignment(BAM::TAlignment& alignment) = 0;

public:
	virtual ~TParsedBamTraverser() = default;
};

enum class ReadType : bool {Filtered, Parsed};

template<ReadType Type>
using TBamReadTraverser = std::conditional_t<Type == ReadType::Parsed, TParsedBamTraverser, void>;
} // namespace GenomeTasks
#endif
