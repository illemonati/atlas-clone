#ifndef TSIMULATORALLELEINDEX_H_
#define TSIMULATORALLELEINDEX_H_

#include <cstddef>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Files/gzstream.h"
#include "genometools/GenotypeTypes.h"

namespace Simulations {
class TSimulatorAlleleIndex {
private:
	size_t nextIndex = 0;
	genometools::Base refBase = genometools::Base::N;
	genometools::Base indexToBase[4];

public:
	coretools::TStrongArray<int, genometools::Base, 4> index;
	coretools::TStrongArray<bool, genometools::Base, 4> used;

	void clear(const genometools::Base &ref) noexcept {
		used.fill(false);
		used[ref]  = true;
		index[ref] = 0;
		nextIndex  = 1;
		refBase    = ref;
	}

	void add(const genometools::Base &base) noexcept {
		if (!used[base]) {
			used[base]             = true;
			index[base]            = nextIndex;
			indexToBase[nextIndex] = base;
			++nextIndex;
		}
	}

	void writeRefAltToVCF(gz::ogzstream &VCF) {
		VCF << refBase << '\t';
		if (nextIndex == 1) // no alt
			VCF << ".";
		else {
			VCF << indexToBase[1];
			for (size_t i = 2; i < nextIndex; ++i) VCF << ',' << indexToBase[i];
		}
	}
};

}

#endif
