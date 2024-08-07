#include "TGLFIndex.h"

#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/TInputFile.h"

namespace GLF {

std::string TGLFIndex::_getIndexFileName(std::string_view FileName){
	std::string name(coretools::str::readBeforeLast(FileName, '.'));
	return name + ".idx";
}

void TGLFIndex::clear(){
	_chrs.clear();
	_posInFile.clear();
}

void TGLFIndex::addChromosme(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PosInFile){
	if(Ploidy != 1 && Ploidy != 2)
		UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Ploidy, ")!");
	_chrs.appendChromosome(Name, Length, Ploidy);
	_posInFile.push_back(PosInFile);
}

void TGLFIndex::addChromosme(const genometools::TChromosome& Chr, uint64_t PosInFile){
	if(Chr.ploidy() != 1 && Chr.ploidy() != 2)
		UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Chr.ploidy(), ")!");
	_chrs.appendChromosome(Chr);
	_posInFile.push_back(PosInFile);
}

void TGLFIndex::readChromosomes(std::string_view GLFFilename) {
	clear();
	for(auto in = coretools::TInputFile(_getIndexFileName(GLFFilename), coretools::FileType::Header); !in.empty(); in.popFront()){
		_chrs.appendChromosome(in.get("Chr"), in.get<uint32_t>("Length"), in.get<uint8_t>("Ploidy"));
		_posInFile.push_back(in.get<uint64_t>("PosInFile"));
	}
}

void TGLFIndex::writeChromosmes(std::string_view GLFFilename){
	coretools::TOutputFile out(_getIndexFileName(GLFFilename), {"Chr", "Length", "Ploidy", "PosInFile"});
	for(auto c: _chrs){
		out.writeln(c.name(), c.length(), c.ploidy(), _posInFile[c.refID()]);
	}
}

void TGLFIndex::checkChromosome(size_t LastRefID, std::string_view Name, uint32_t Length, uint8_t Ploidy) {
	size_t oldRefID = LastRefID;

	// forward
	while (LastRefID < _chrs.size() && _chrs[LastRefID].name() != Name) {++LastRefID;}
	if (LastRefID == _chrs.size()) {
		LastRefID = 0;
		// wrap around
		while (LastRefID < oldRefID && _chrs[LastRefID].name() != Name) {++LastRefID;}
	}
	if (_chrs[LastRefID].name() !=  Name)
		UERROR("Chromosome '", Name, "' not found in GLF index!");

	if(_chrs[LastRefID].length() != Length)
		UERROR("Length of chromosome '", Name, "' does not match length given in GLF index!");

	if(_chrs[LastRefID].ploidy() != Ploidy)
		UERROR("Ploidy of chromosome '", Name, "' does not match ploidy given in GLF index!");
}

bool TGLFIndex::hasSameChromosomes(const TGLFIndex &Other) const {
	// checks if two TGlfIndexFiles contain the same chromosomes in terms of order, names and lengths
	// ploidy may be different (e.g. X for males and females), as well as posInFile, which depends on the amount of data
	if(_chrs.size() != Other._chrs.size()) return false;

	auto a = _chrs.begin();
	auto b = Other._chrs.begin();
	while(a != _chrs.end()){
		if(a->name() != b->name() || a->length() != b->length()) return false;
		++a; ++b;
	}
	return true;
}
}
