/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include <array>
#include <map>
#include <string>
#include <vector>
#include <zlib.h>

#include "genometools/GenotypeTypes.h"
#include "genometools/PhredProbabilityTypes.h"
#include "coretools/Containers/TDualStrongArray.h"
#include "GenotypeData.h"
#include "genometools/VCF/TVcfWriter.h"
#include "genometools/GenomePositions/TChromosomes.h"

namespace GLF {
using genometools::Ploidy;
using TGLFLikelihoods = coretools::TDualStrongArray<genometools::HighPrecisionPhredIntProbability, genometools::Base, genometools::Genotype, 4, 10, Ploidy>;

//----------------------------------------------------
// TGlfChromosome
// struct to store info on chromosome
// TODO: derive from TChromosome??
//----------------------------------------------------
class TGlfChromosome {
private:
	std::string _name = "";
	uint32_t _refId = 0;
	uint32_t _length = 0;
	bool _isHaploid = false;

public:
	TGlfChromosome() = default;
	TGlfChromosome(const std::string &Name, uint32_t Length, uint8_t Ploidy);

	std::string name() const { return _name; };
	uint32_t refId() const noexcept { return _refId; };
	uint32_t length() const noexcept { return _length; };
	bool isHaploid() const noexcept { return _isHaploid; };
	uint8_t ploidy() const noexcept { return 2 - _isHaploid; };
	size_t numLikelihoodValues() const noexcept {
		std::array<size_t, 2> N{10, 4};
		return N[_isHaploid];
	};

	void update(std::string_view Name, uint16_t RefId, uint32_t Length, uint8_t Ploidy);
	void clear();
};

//----------------------------------------------------
// TGlfHandle
//----------------------------------------------------
class TGlfHandle {
protected:
	std::string _filename;
	gzFile _gzfp         = nullptr;
	std::string _version = "GLF3";
	uint64_t _positionInFile = 0;
	TGlfChromosome _curChr;

public:
	TGlfHandle()                              = default;
	TGlfHandle(const TGlfHandle &)            = delete;
	TGlfHandle &operator=(const TGlfHandle &) = delete;
	TGlfHandle(TGlfHandle &&)                 = default;
	TGlfHandle &operator=(TGlfHandle &&)      = default;

	virtual ~TGlfHandle() {close();}

	void close() {
		if (_gzfp) {
			gzclose(_gzfp);
			_gzfp = nullptr;
		}
	};

	std::string name() const { return _filename; };

	std::string chr() const { return _curChr.name(); };

	uint32_t refId() const noexcept { return _curChr.refId(); };

	uint32_t chrLength() const noexcept { return _curChr.length(); };

	bool chrIsHaploid() const noexcept { return _curChr.isHaploid(); };

	uint8_t chrNumLikelihoodValues() const noexcept { return _curChr.numLikelihoodValues(); };
};

//----------------------------------------------------
// TGlfWriter
//----------------------------------------------------
class TGlfWriter : public TGlfHandle {
private:
	long _oldPos         = 0;
	std::string _header;

	void _writeHeader();

	template<typename T> void _write(T var) { _positionInFile += gzwrite(_gzfp, &var, sizeof(T)); };
	void _write(const void *buf, size_t len) { _positionInFile += gzwrite(_gzfp, buf, len); };

public:
	TGlfWriter() = default;
	TGlfWriter(const std::string &Filename, const genometools::TChromosomes &Chrs) {
		open(Filename, Chrs, "");
	};

	// open & close streams
	void open(const std::string &Filename, const genometools::TChromosomes &Chrs, const std::string &Header = "");
	void newChromosome(const genometools::TChromosome &chromosome);
	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual,
		       const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods);
};

//----------------------------------------------------
// TGlfReader
//----------------------------------------------------
class TGlfReader : public TGlfHandle {
private:
	// file parsing
	std::string _glfFileVersion;
	bool _eof = true;

	// about site
	int _recordType = 99;
	uint32_t _position = 0;
	uint16_t _depth = 0;
	int _RMS_mappingQual = 0;
	TGLFLikelihoods _genotypeLikelihoodsGLF;

	// about chromosomes
	std::map<uint32_t, TGlfChromosome> _chromosomesAlreadyParsed;
	genometools::TChromosomes _chromosomes;

	// initializing
	void _open();

	// read
	template<typename T> int _read(T *buf, size_t len) {
		const int _lenRead = gzread(_gzfp, buf, len);
		return _lenRead;
	};

	template<typename T> T _read() {
		T t;
		const int _lenRead = gzread(_gzfp, &t, sizeof(T));
		if (!_lenRead) UERROR("could not read file ", name(), "!");
		_positionInFile += _lenRead;
		return t;
	};

	bool _readChr(bool storePrevious = true);
	bool _readRecordType();
	void _readSNPRecord();
	bool _jumpToEndOfChr();

public:
	TGlfReader() = default;
	TGlfReader(const std::string &Filename) {
		open(Filename);
	};

	// get details
	bool eof() const { return _eof; };
	const genometools::TChromosomes& chromosomes();
	TGlfChromosome *pointerToChr(uint32_t refId);	
	uint32_t position() const { return _position; };
	uint16_t depth() const { return _depth; };
	const TGLFLikelihoods &genotypeLikelihoodsGLF() const { return _genotypeLikelihoodsGLF; };	
	TGlfChromosome *curChr() {return &_curChr;};
	
	// open file and parse header
	void setFilename(const std::string &Filename);
	void open(const std::string &Filename);
	void rewind();

	// parse file
	bool readNext();
	bool jumpToNextChr();
	bool readNextWindow(std::vector<TGLFLikelihoods> &genoLikelihoods, uint32_t refId, uint32_t start, uint32_t end);

	// printing
	void printChr();
	void printSite();
	void printToEnd();
};

struct TGLFPrinter {
	void run();
};

}; // namespace GLF

#endif /* TGLF_H_ */
