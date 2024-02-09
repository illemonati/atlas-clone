
#ifndef TSIMULATORMUTATIONTABLE_H_
#define TSIMULATORMUTATIONTABLE_H_

#include "GenotypeData.h"
#include "coretools/Containers/TStrongArray.h"
#include "genometools/GenotypeTypes.h"

namespace Simulations {
class TSimulatorMutationtable {
private:
	coretools::TStrongArray<coretools::TStrongArray<double, genometools::Base>, genometools::Base> _mutTable;
	void _normalizeAndMakeCumulative();
public:
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq);
	TSimulatorMutationtable(const GenotypeLikelihoods::TBaseProbabilities &baseFreq, double theta);

	const coretools::TStrongArray<double, genometools::Base> &operator[](genometools::Base base) const noexcept { return _mutTable[base]; }
};
}

#endif
