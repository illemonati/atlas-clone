#ifndef TGLFWRITER_H_
#define TGLFWRITER_H_

#include "GenotypeData.h"
#include "TGLFIndex.h"
#include "coretools/Files/TWriter.h"
#include <string>

namespace GLF {

class TGLFWriter {
private:
	long _oldPos         = 0;
	std::string _header;

	TGLFIndex _index;
	std::unique_ptr<coretools::TWriter> _writer;

	void _writeHeader(std::string_view Header);

public:
	TGLFWriter() = default;
	TGLFWriter(std::string_view Filename, const genometools::TChromosomes &Chrs) {
		open(Filename, Chrs, "");
	};
	~TGLFWriter(){ close(); };

	// open & close streams
	void open(std::string_view Filename, const genometools::TChromosomes &Chrs, std::string_view Header = "");
	void newChromosome(const genometools::TChromosome &chromosome);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual,
			   const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods);
	void close();
};

}

#endif
