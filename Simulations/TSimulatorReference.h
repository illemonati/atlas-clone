
#ifndef TSIMULATORREFERENCE_H_
#define TSIMULATORREFERENCE_H_

#include <fstream>
#include <vector>

#include "genometools/Genotypes/Base.h"
#include "genometools/TFastaWriter.h"

namespace Simulations {

class TSimulatorReference {
private:
	genometools::TFastaWriter _writer;
	std::vector<genometools::Base> _ref;

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
