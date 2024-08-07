#ifndef TGLFINDEX_H_
#define TGLFINDEX_H_

#include "genometools/GenomePositions/TChromosomes.h"
 
namespace GLF {
class TGLFIndex{
private:
	genometools::TChromosomes _chrs;
	std::vector<uint64_t> _posInFile;

	std::string _getIndexFileName(std::string_view FileName);

public:
	void clear();

	void addChromosme(std::string_view Name, uint32_t Length, uint8_t Ploidy, uint64_t PosInFile);
	void addChromosme(const genometools::TChromosome& Chr, uint64_t PosInFile);

	void writeChromosmes(std::string_view GLFFilename);
	void readChromosomes(std::string_view GLFFilename);

	void checkChromosome(size_t LastRefID, std::string_view Name, uint32_t Length, uint8_t Ploidy);

	// compare
	bool hasSameChromosomes(const TGLFIndex& Other) const;

	// getters do not check if chromosomes were initialized!
	const genometools::TChromosomes& chromosomes() const noexcept { return _chrs; };
	size_t lastChrNumLikelihoodValues() const noexcept {
		std::array<size_t, 3> N{0, 4, 10}; // for ploidy 0, 1 and 2
		return N[_chrs.back().ploidy()];
	}

	size_t chrNumLikelihoodValues(size_t refID) const noexcept {
		std::array<size_t, 3> N{0, 4, 10}; // for ploidy 0, 1 and 2
		return N[_chrs[refID].ploidy()];
	}

	uint64_t positionInFile(size_t RefID) const noexcept { return _posInFile[RefID]; }
	size_t size() const noexcept {return _chrs.size();}

};
}

#endif
