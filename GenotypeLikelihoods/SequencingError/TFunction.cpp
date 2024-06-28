#include "TFunction.h"
#include "coretools/Main/TLog.h"

namespace GenotypeLikelihoods::SequencingError {

void TFunction::log() const {
	using coretools::instances::logfile;
	using coretools::str::toString;
	logfile().write();
	constexpr size_t Nmax = 3;
	const size_t Nbeta    = std::distance(begin(), end());
	if (Nbeta == 0) {
		logfile().list(typeString(), ": []");
		return;
	}

	std::string out = "[";
	if (Nbeta <= Nmax * 2) {
		for (auto p : *this) out.append(toString(p, ", "));
	} else {
		for (size_t i = 0; i < Nmax; ++i) out.append(toString(*(begin() + i), ", "));
		out.append("...,");
		const auto pStart = end() - Nmax;
		for (size_t i = 0; i < Nmax; ++i) out.append(toString(*pStart + i, ", "));
	}
	out.pop_back();
	out.back() = ']';
	logfile().list(typeString(), ": ", out);
}
} // namespace GenotypeLikelihoods::SequencingError
