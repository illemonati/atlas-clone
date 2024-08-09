#ifndef TGLFREADER_H_
#define TGLFREADER_H_

#include "TGLF.h"
#include "TGLFIndex.h"
#include "coretools/Files/TReader.h"
#include "genometools/GenomePositions/TGenomePosition.h"

namespace GLF {

class TGLFReader {
private:
	// file parsing
	bool _eof      = true;
	bool _hasIndex = false;

	// about site
	TGLFLikelihoods _genotypeLikelihoodsGLF;
	genometools::TGenomePosition _position;
	int _recordType          = 99;
	uint16_t _depth          = 0;
	uint8_t _RMS_mappingQual = 0;

	TGLFIndex _index;
	std::unique_ptr<coretools::TReader> _reader = std::make_unique<coretools::TNoReader>();

	bool _readChr();
	bool _readRecordType();
	void _readSNPRecord();

public:
	TGLFReader() = default;
	TGLFReader(const std::string &Filename, bool HasIndex = true) {
		open(Filename, HasIndex);
	}

	bool empty() const noexcept { return _eof; };
	void popFront();
	const TGLFLikelihoods &front() const noexcept { return _genotypeLikelihoodsGLF; };

	// get details
	const genometools::TChromosomes& chromosomes() const;
	const genometools::TChromosome& curChromosome() const;
	size_t lastRefID() const;
	uint32_t refID() const noexcept { return _position.refID(); };
	const TGLFIndex& index() const noexcept { return _index; };
	genometools::TGenomePosition position() const noexcept { return _position; };
	uint16_t depth() const noexcept { return _depth; };

	const std::string& name() const {return _reader->name();}

	// open file and parse header
	void open(const std::string &Filename, bool HasIndex = true);
	void rewind();

	// parse file
	bool jumpToChr(uint32_t RefID, bool OnlyForward=true);
	bool jumpToNextChr();
	bool jumpToPositionOrBeyond(const genometools::TGenomePosition& Position, bool OnlyForward=true);
	bool readNextWindow(std::vector<TGLFLikelihoods> &GenoLikelihoods, genometools::TGenomeWindow Window);

	// printing
	void printChr();
	void printSite();
	void printToEnd();
	void writeIndex();

	const genometools::TChromosome& curChr() const noexcept {return _index.chromosomes()[refID()];}
};

}
#endif
