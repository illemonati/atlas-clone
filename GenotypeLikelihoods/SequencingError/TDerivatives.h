#ifndef GENOTYPELIKELIHOODS_TDERIVATIVES_H_
#define GENOTYPELIKELIHOODS_TDERIVATIVES_H_

#include <cstddef>

namespace GenotypeLikelihoods {
namespace SequencingError {

struct T1stDerivative {
	size_t index;
	double derivative;
	T1stDerivative(size_t Index, double Derivative) : index(Index), derivative(Derivative) {}
};

struct T2ndDerivative {
	size_t index1;
	size_t index2;
	double derivative;
	T2ndDerivative(size_t Index1, size_t Index2, double Derivative) : index1(Index1), index2(Index2), derivative(Derivative) {}
};
} // namespace SequencingError
} // namespace GenotypeLikelihoods
#endif
