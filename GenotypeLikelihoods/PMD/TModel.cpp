#include "TModel.h"
#include "TPerReadGroup.h"
#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "coretools/Types/probability.h"
#include "genometools/GenotypeTypes.h"
#include <utility>

namespace GenotypeLikelihoods::PMD {
using genometools::Base;
using namespace coretools::str;

TBaseProbabilities TNoPMD::massFunction(genometools::Genotype g, const TBaseLikelihoods &baseLikelihoodsNoPMD) {
	using namespace genometools;
	const Base a = first(g);
	const Base b = second(g);
	TBaseLikelihoods mf{0};
	mf[a] = baseLikelihoodsNoPMD[a];
	mf[b] = baseLikelihoodsNoPMD[b];
	return TBaseProbabilities::normalize(mf);
}


TModel *makeType(std::string_view pmdString) {
	// split by ':'
	std::vector<std::string> details;
	fillContainerFromString(pmdString, details, ":");

	// switch type
	if (details[0] == TNoPMD::name) return new TNoPMD;
	if (details[0] == TPerReadGroup<Strand::Single>::name) return new TPerReadGroup<Strand::Single>(details);
	if (details[0] == TPerReadGroup<Strand::Double>::name) return new TPerReadGroup<Strand::Double>(details);

	UERROR("Cannot initialize PMD: unknown PMD type '", details[0], "'!\nUse ", TNoPMD::name, " or ",
		   TPerReadGroup<Strand::Single>::name, " or ", TPerReadGroup<Strand::Double>::name, ".");
}
} // namespace GenotypeLikelihoods
