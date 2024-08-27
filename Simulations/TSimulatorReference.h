
#ifndef TSIMULATORREFERENCE_H_
#define TSIMULATORREFERENCE_H_

#include <fstream>
#include <vector>

#include "genometools/Genotypes/Base.h"

namespace Simulations {

class TSimulatorReference {
private:

	// fasta file
	std::ofstream _fasta;
	std::ofstream _fastaIndex;
	long _oldOffset = 0;
	std::string _filename;

	// reference storage
	std::vector<genometools::Base> _ref;
	long _chrLength      = 0;
	std::string _chrName = "";
	bool _needsWriting   = false;

	void _allocateStorage(long length);
	void _closeFastaFile();
	void _writeRefToFasta();

public:
	TSimulatorReference() = default;
	TSimulatorReference(std::string Filename);
	~TSimulatorReference();
	void open(std::string Filename);
	void setChr(std::string ChrName, long ChrLength);
	//	void simulateReferenceSequenceCurChromosome(float* cumulBaseFreq);

	genometools::Base &operator[](size_t index) { return _ref[index]; }
	const genometools::Base &operator[](size_t index) const { return _ref[index]; }
	const std::vector<genometools::Base> & reference() const {return _ref;}
};
}

#endif
