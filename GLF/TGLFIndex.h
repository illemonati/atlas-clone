#ifndef TGLFINDEX_H_
#define TGLFINDEX_H_

#include "genometools/GenomePositions/TChromosomes.h"
#include <cstdint>
 
namespace GLF {
class TGLFIndex{
private:
	struct TEntry {
		std::string name;
		size_t length;
		uint8_t ploidy;
		size_t position;
		TEntry(std::string_view Name, size_t Length, uint8_t Ploidy, size_t Position)
			: name(Name), length(Length), ploidy(Ploidy), position(Position) {};
	};
	genometools::TChromosomes _chrs;
	std::vector<TEntry> _entries;

	std::string _getIndexFileName(std::string_view FileName);

public:
	void clear();

	void addChromosme(std::string_view Name, size_t Length, uint8_t Ploidy, size_t PosInFile);
	void addChromosme(const genometools::TChromosome& Chr, uint64_t PosInFile);

	void writeChromosmes(std::string_view GLFFilename);
	void readChromosomes(std::string_view GLFFilename);

	void checkChromosome(size_t RefID, std::string_view Name, uint32_t Length, uint8_t Ploidy);

	// compare
	bool hasSameChromosomes(const TGLFIndex& Other) const;

	// getters do not check if chromosomes were initialized!
	const genometools::TChromosomes& chromosomes() const noexcept { return _chrs; };
	size_t lastChrNumLikelihoodValues() const noexcept {
		assert(size() > 0);
		std::array<size_t, 3> N{0, 4, 10}; // for ploidy 0, 1 and 2
		return N[_entries.back().ploidy];
	}

	size_t chrNumLikelihoodValues(size_t RefID) const noexcept {
		assert(RefID < size());
		std::array<size_t, 3> N{0, 4, 10}; // for ploidy 0, 1 and 2
		return N[_entries[RefID].ploidy];
	}

	const std::string& name(size_t RefID) const noexcept {
		assert(RefID < size());
		return _entries[RefID].name;
	}

	uint64_t positionInFile(size_t RefID) const noexcept {
		assert(RefID < size());
		return _entries[RefID].position;
	}
	uint64_t length(size_t RefID) const noexcept {
		assert(RefID < size());
		return _entries[RefID].length;
	}
	uint8_t ploidy(size_t RefID) const noexcept {
		assert(RefID < size());
		return _entries[RefID].ploidy;
	}
	size_t size() const noexcept {return _chrs.size();}

	bool empty() const noexcept {return size() == 0;}

};
}

#endif
