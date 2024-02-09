#ifndef SEQUENCINGERROR_TINTERCEPT_H_
#define SEQUENCINGERROR_TINTERCEPT_H_

#include "SequencingError/TDerivatives.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods::SequencingError {

class TIntercept final  {
private:
	double _beta    = 0.;

public:
	static constexpr std::string_view name = "intercept";

	constexpr size_t numParameters() const noexcept { return 1; }

	constexpr double &beta() noexcept {return _beta;}
	constexpr double beta() const noexcept {return _beta;}

	constexpr double getEta() const noexcept { return _beta; }

	double getEta(std::vector<T1stDerivative> &der1) const noexcept {
		der1.emplace_back(0, 1.);
		return _beta;
	}

	std::string typeString() const noexcept { return std::string(name); }

	void addInfo(BAM::RGInfo::TInfo& info) const;
	void log() const;
};
}

#endif
