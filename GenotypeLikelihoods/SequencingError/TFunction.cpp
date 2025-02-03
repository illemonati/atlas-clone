#include "TFunction.h"
#include "coretools/Main/TLog.h"

namespace GenotypeLikelihoods::SequencingError {

void TFunction::log() const {
	using coretools::instances::logfile;
	using coretools::str::toString;
	logfile().write();
	const size_t Nbeta    = std::distance(begin(), end());
	if (Nbeta == 0) {
		logfile().list(typeString(), ": []");
		return;
	}

	std::string out = "[";
	for (auto p : *this) out.append(toString(p, ", "));
	out.pop_back();
	out.back() = ']';
	logfile().list(typeString(), ": ", out);
}
} // namespace GenotypeLikelihoods::SequencingError
