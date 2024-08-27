
#ifndef TSIMULATORMUTATIONTABLE_H_
#define TSIMULATORMUTATIONTABLE_H_

#include "genometools/Genotypes/Containers.h"
#include "coretools/Containers/TStrongArray.h"

namespace Simulations {
class TSimulatorMutationtable {
private:
	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _mutTable;
	void _normalizeAndMakeCumulative();
public:
	TSimulatorMutationtable(const genometools::TBaseProbabilities &baseFreq);
	TSimulatorMutationtable(const genometools::TBaseProbabilities &baseFreq, double theta);

	const coretools::TStrongArray<double, genometools::Base> &operator[](genometools::Base base) const noexcept { return _mutTable[base]; }
};
}

#endif
