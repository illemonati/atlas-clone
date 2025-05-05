
#ifndef TSIMULATORREFERENCE_H_
#define TSIMULATORREFERENCE_H_

#include <fstream>
#include <vector>

#include "genometools/Genotypes/Base.h"

namespace Simulations {

class TSimulatorReference {
private:
	std::ofstream _fasta;
	std::ofstream _fastaIndex;
	long _oldOffset = 0;
	std::string _filename;

	// reference storage
	std::vector<genometools::Base> _ref;
	std::string _chrName = "";

	void _writeRefToFasta();

public:
	TSimulatorReference(std::string_view Filename);
	~TSimulatorReference();

	void setChr(std::string_view ChrName, long ChrLength);

	genometools::Base &operator[](size_t index) { return _ref[index]; }
	const genometools::Base &operator[](size_t index) const { return _ref[index]; }
};
}

#endif
