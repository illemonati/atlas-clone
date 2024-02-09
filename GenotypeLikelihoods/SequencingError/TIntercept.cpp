#include "TIntercept.h"
#include "coretools/Main/TLog.h"
namespace GenotypeLikelihoods::SequencingError {

	void TIntercept::addInfo(BAM::RGInfo::TInfo& info) const {
		 info[name] = _beta;
	}

	void TIntercept::log() const {
		coretools::instances::logfile().list(typeString(), ": ", _beta);
	}
}
