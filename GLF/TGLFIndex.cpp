#include "TGLFIndex.h"

#include "coretools/Files/TInputFile.h"
#include "coretools/Files/TOutputFile.h"

namespace GLF {

std::string TGLFIndex::_getIndexFileName(std::string_view FileName) {
	std::string name(coretools::str::readBeforeLast(FileName, '.'));
	return name + ".idx";
}

void TGLFIndex::clear() {
	_chrs.clear();
	_entries.clear();
}

void TGLFIndex::addChromosme(std::string_view Name, size_t Length, uint8_t Ploidy, size_t PosInFile) {
	if (Ploidy != 1 && Ploidy != 2) UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Ploidy, ")!");
	_chrs.appendChromosome(Name, Length, Ploidy);
	_entries.emplace_back(Name, Length, Ploidy, PosInFile);
}

void TGLFIndex::addChromosme(const genometools::TChromosome &Chr, uint64_t PosInFile) {
	if (Chr.ploidy() != 1 && Chr.ploidy() != 2)
		UERROR("Currently GLFs only support ploidies 1 and 2 (not ", Chr.ploidy(), ")!");
	_chrs.appendChromosome(Chr);
	_entries.emplace_back(Chr.name(), Chr.length(), Chr.ploidy(), PosInFile);
}

void TGLFIndex::readChromosomes(std::string_view GLFFilename) {
	clear();
	for (auto in = coretools::TInputFile(_getIndexFileName(GLFFilename), coretools::FileType::Header); !in.empty();
	     in.popFront()) {
		_chrs.appendChromosome(in.get("Chr"), in.get<size_t>("Length"), in.get<uint8_t>("Ploidy"));
		_entries.emplace_back(in.get("Chr"), in.get<size_t>("Length"), in.get<uint8_t>("Ploidy"),
		                      in.get<size_t>("PosInFile"));
	}
}

void TGLFIndex::writeChromosmes(std::string_view GLFFilename) {
	coretools::TOutputFile out(_getIndexFileName(GLFFilename), {"Chr", "Length", "Ploidy", "PosInFile"});
	for (const auto &e : _entries) { out.writeln(e.name, e.length, e.ploidy, e.position); }
}

void TGLFIndex::checkChromosome(size_t RefID, std::string_view Name, uint32_t Length, uint8_t Ploidy) {
	if (RefID >= _chrs.size()) {
		UERROR("Chromosome ", RefID, ", named ", Name, " in GLF file is a higher number than the ", _chrs.size(), " chromosomes in index file!");
	}

	if (_chrs[RefID].name() != Name) {
		UERROR("Chromosome ", RefID, " is named ", Name, " in GLF file but ", _chrs[RefID].name(), " in index file!");
	}
	if (_chrs[RefID].length() != Length) {
		UERROR("Chromosome ", RefID, " has length ", Length, " in GLF file but ", _chrs[RefID].length(),
		       " in index file!");
	}
	if (_chrs[RefID].ploidy() != Ploidy) {
		UERROR("Chromosome ", RefID, " has ploidy ", Ploidy, " in GLF file but ", _chrs[RefID].ploidy(),
		       " in index file!");
	}
}

bool TGLFIndex::hasSameChromosomes(const TGLFIndex &Other) const {
	// checks if two TGlfIndexFiles contain the same chromosomes in terms of order, names and lengths
	// ploidy may be different (e.g. X for males and females), as well as posInFile, which depends on the amount of data
	if (size() != Other.size()) return false;

	for (size_t i = 0; i < size(); ++i) {
		const auto &a = chromosomes()[i];
		const auto &b = Other.chromosomes()[i];
		if (a.name() != b.name() || a.length() != b.length()) return false;
	}
	return true;
}
} // namespace GLF
