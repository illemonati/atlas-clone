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

//--------------------------------------------
// TGlfIndexFile
//--------------------------------------------
class TGlfIndexFile{
private:
	std::string _indexFilename;
	coretools::TOutputFile _indexFileOutput;

	genometools::TChromosomes _chrs;
	std::vector<uint64_t> _posInFile;
	size_t _curChr;	

public:
	TGlfIndexFile(){};
	~TGlfIndexFile() = default;
	
	void clear();

	void addChromosme(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PosInFile);
	void addChromosme(const genometools::TChromosome& Chr, uint64_t PosInFile);

	void writeChromosmes(std::string_view Filename);
	void readChromosomes(std::string_view Filename);

	void jumpToNextChromosome();
	void seekChr(uint32_t RefID);
	void setCurChromosomeToAndCheck(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PositionInFile);

	// compare
	bool hasSameChromosomes(const TGlfIndexFile& Other) const;

	// getters do not check if chromosomes were initialized!
	const genometools::TChromosomes& chromosomes(){ return _chrs; };
	const genometools::TChromosome& curChr(){ return _chrs[_curChr]; };
	size_t curChrNumLikelihoodValues() const noexcept {
		std::array<size_t, 2> N{10, 4};
		return N[_chrs[_curChr].ploidy()];
	};	
	uint64_t curChrPositionInFile(){ return _posInFile[_curChr]; };
	uint64_t nextChrPositionInFile();
	bool curChrIsLast(){ return _curChr == _chrs.size() - 1; };
};

//----------------------------------------------------
// TGlfHandle
//----------------------------------------------------
class TGlfHandle {
protected:
	std::string _filename;
	gzFile _gzfp         = nullptr;
	
	std::string _version = "GLF2";
	uint64_t _positionInFile = 0;

	TGlfIndexFile _index;

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
	//uint8_t chrNumLikelihoodValues() const noexcept { return _curChr.numLikelihoodValues(); };
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
	bool _hasIndex = false;

	// about site
	int _recordType = 99;
	uint32_t _position = 0;
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
	TGlfReader() = default;
	TGlfReader(const std::string &Filename) {
		open(Filename);
	};

	// get details
	bool eof() const { return _eof; };
	const genometools::TChromosomes& chromosomes();
	const genometools::TChromosome& curChromosome();	
	uint32_t refId(){ return _index.curChr().refID(); };
	const TGlfIndexFile& index(){ return _index; };	
	uint32_t position() const { return _position; };
	uint16_t depth() const { return _depth; };
	const TGLFLikelihoods &genotypeLikelihoodsGLF() const { return _genotypeLikelihoodsGLF; };		
	
	// open file and parse header
	void setFilename(const std::string &Filename);
	void open(const std::string &Filename, bool HasIndex = true);
	void rewind();

	// parse file
	bool readNext();
	bool jumpToChr(uint32_t RefID);
	bool jumpToNextChr();
	bool readNextWindow(std::vector<TGLFLikelihoods> &genoLikelihoods, uint32_t refId, uint32_t start, uint32_t end);

	// printing
	void printChr();
	void printSite();
	void printToEnd();
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
