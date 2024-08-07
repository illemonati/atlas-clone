/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include <string>
#include <vector>
#include <zlib.h>

#include "TGLFIndex.h"

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Containers/TDualStrongArray.h"
#include "genometools/VCF/TVcfWriter.h"
#include "genometools/GenomePositions/TChromosomes.h"

namespace GLF {
using genometools::Ploidy;
using TGLFLikelihoods = coretools::TDualStrongArray<genometools::HighPrecisionPhredIntProbability, genometools::Base, genometools::Genotype, 4, 10, Ploidy>;

constexpr std::string_view version() noexcept { return "GLF2"; }

//----------------------------------------------------
// TGlfHandle
//----------------------------------------------------
class TGlfHandle {
protected:
	std::string _filename;
	gzFile _gzfp         = nullptr;

	std::string _version = "GLF2";
	uint64_t _positionInFile = 0;

	TGLFIndex _index;

public:
	TGlfHandle()                              = default;
	TGlfHandle(const TGlfHandle &)            = delete;
	TGlfHandle &operator=(const TGlfHandle &) = delete;
	TGlfHandle(TGlfHandle &&)                 = default;
	TGlfHandle &operator=(TGlfHandle &&)      = default;

	virtual ~TGlfHandle() {close();}

	virtual void close() {
		if (_gzfp) {
			gzclose(_gzfp);
			_gzfp = nullptr;
		}
	};

	std::string name() const { return _filename; };
	//uint8_t chrNumLikelihoodValues() const noexcept { return _curChr.numLikelihoodValues(); };
};

//----------------------------------------------------
// TGlfReader
//----------------------------------------------------
class TGLFReader : public TGlfHandle {
private:
	// file parsing
	std::string _glfFileVersion;
	bool _eof = true;
	bool _hasIndex = false;

	// about site
	int _recordType = 99;
	genometools::TGenomePosition _position;
	uint16_t _depth = 0;
	int _RMS_mappingQual = 0;
	TGLFLikelihoods _genotypeLikelihoodsGLF;

	// initializing
	void _open(bool HasIndex = true);

	// read
	template<typename T> int _read(T *buf, size_t len) {
		const int _lenRead = gzread(_gzfp, buf, len);
		_positionInFile += _lenRead;
		return _lenRead;
	};

	template<typename T> T _read() {
		T t;
		const int _lenRead = gzread(_gzfp, &t, sizeof(T));
		if (!_lenRead) UERROR("could not read file ", name(), "!");
		_positionInFile += _lenRead;
		return t;
	};

	bool _readChr();
	bool _readRecordType();
	void _readSNPRecord();

public:
	TGLFReader() = default;
	TGLFReader(const std::string &Filename, bool HasIndex = true) {
		open(Filename, HasIndex);
	};

	// get details
	bool eof() const { return _eof; };
	const genometools::TChromosomes& chromosomes();
	const genometools::TChromosome& curChromosome();
	size_t lastRefID();
	uint32_t refID() const noexcept { return _position.refID(); };
	const TGLFIndex& index() const noexcept { return _index; };
	genometools::TGenomePosition position() const noexcept { return _position; };
	uint16_t depth() const noexcept { return _depth; };
	const TGLFLikelihoods &genotypeLikelihoodsGLF() const noexcept { return _genotypeLikelihoodsGLF; };

	// open file and parse header
	void setFilename(const std::string &Filename);
	void open(const std::string &Filename, bool HasIndex = true);
	void rewind();

	// parse file
	bool readNext();
	bool jumpToChr(uint32_t RefID, bool OnlyForward=true);
	bool jumpForwardToChr(uint32_t RefID);
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


//--------------------------------------------
// Tasks
//--------------------------------------------
struct TGLFPrinter {
	void run();
};

struct TGLFIndexer {
	void run();
};

}; // namespace GLF

#endif /* TGLF_H_ */
