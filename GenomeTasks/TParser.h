#ifndef TPARSER_H_
#define TPARSER_H_

#include "genometools/TFastaReader.h"

#include "TBaseFilter.h"
#include "TGenome.h"
#include "TErrorModels.h"

namespace GenomeTasks {

class TParser {
	bool _trimReads;
	int _trim3;
	int _trim5;

	// filters
	TQualityFilter _qualityFilter;
	TContextFilter _contextFilter;

	genometools::TFastaReader _reference;
	GenotypeLikelihoods::TErrorModels _errorModels;
public:
	TParser(TGenome& genome);

	void apply(BAM::TAlignment &alignment);

	void openReference(bool required = false);
	const genometools::TFastaReader& reference() const noexcept {return _reference;};
	const GenotypeLikelihoods::TErrorModels& errorModels() const noexcept {return _errorModels;};
};
}

#endif
