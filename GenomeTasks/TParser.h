#ifndef TPARSER_H_
#define TPARSER_H_

#include "genometools/TFastaReader.h"

#include "TBaseFilter.h"
#include "TGenome.h"

namespace GenomeTasks {

class TParser {
	bool _trimReads;
	int _trim3;
	int _trim5;

	TQualityFilter _qualityFilter;
	TContextFilter _contextFilter;

	genometools::TFastaReader _reference;

public:
	TParser();

	void fill(const TGenome& genome, BAM::TAlignment& alignment) const;
	void openReference(bool required = false);
	const genometools::TFastaReader& reference() const noexcept {return _reference;};
};
}

#endif
