/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include <array>
#include <cstdint>
#include <iomanip>
#include <memory>
#include <stddef.h>
#include <string>
#include <vector>

#include "PhredProbabilityTypes.h"
#include "devtools.h"
#include "stringFunctions.h"
#include "mathFunctions.h"
#include "probability.h"

#include <armadillo>

namespace GenotypeLikelihoods {
namespace SequencingError {

struct T1stDerivative{
       size_t index;
       double derivative;
};

struct T2ndDerivative{
       size_t index1;
       size_t index2;
       double derivative;
};
//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TFunction {
private:
	size_t _firstParameterIndex;

protected:
	virtual double *_begin() noexcept              = 0;
	virtual double *_end() noexcept                = 0;
	virtual const double *_cbegin() const noexcept = 0;
	virtual const double *_cend() const noexcept   = 0;
	virtual double *_obegin() noexcept             = 0;
	virtual double *_oend() noexcept               = 0;

	void _initializeValues(const std::vector<std::string> &betas);

public:
	TFunction(uint16_t FirstParameterIndex = 0) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction() = default;

	// non-virtuals
	void proposeNewParameters(const arma::mat &JxF, uint16_t &index, double lambda) noexcept;
	void rejectProposedParameters() noexcept;
	constexpr size_t firstParameterIndex() const noexcept { return _firstParameterIndex; };

	// virtuals
	virtual uint16_t numParameters() const noexcept               = 0;
	virtual uint16_t numNonZeroFirstDerivatives() const noexcept  = 0;
	virtual uint16_t numNonZeroSecondDerivatives() const noexcept = 0;

	virtual T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept           = 0;
	virtual T2ndDerivative get2ndDerivatives(uint16_t val, size_t i, size_t j) const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual bool checkOrInitValueRange(const std::vector<uint16_t> &values) = 0;

	// estimation
	virtual double getEtaTerm(uint16_t val) const noexcept                = 0;
	virtual double adjustParametersPostEstimation() noexcept              = 0;
	virtual std::string typeString() const noexcept                       = 0;
	std::string modelString() const;
};

//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
class TIntercept : public TFunction {
private:
	double _beta    = 0.;
	double _oldBeta = 0.;

protected:
	double *_begin() noexcept override { return &_beta; }
	double *_end() noexcept override { return &_beta + 1; }
	const double *_cbegin() const noexcept override { return &_beta; }
	const double *_cend() const noexcept override { return &_beta + 1; }
	double *_obegin() noexcept override { return &_oldBeta; }
	double *_oend() noexcept override { return &_oldBeta + 1; }

public:
	static inline const std::string name = "intercept";

	TIntercept(uint16_t FirstParameterIndex, const std::string &beta)
		: TFunction(FirstParameterIndex), _beta(coretools::str::convertStringCheck<double>(beta)){};

	uint16_t numParameters() const noexcept override { return 1; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	bool checkOrInitValueRange(const std::vector<uint16_t> &) noexcept override { return true; };

	constexpr double intercept() const noexcept { return _beta; }
	constexpr double &intercept() noexcept { return _beta; }

	double getEtaTerm(uint16_t = 0) const noexcept override { return _beta; };
	T1stDerivative get1stDerivatives() const noexcept { return {firstParameterIndex(), 1.}; };
	T1stDerivative get1stDerivatives(uint16_t, size_t) const noexcept override { return {firstParameterIndex(), 1.}; };
	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};
	double adjustParametersPostEstimation() noexcept override { return 0.; }

	virtual std::string typeString() const noexcept override { return name; };
};

namespace impl {
struct TNoTransform {
	static inline const std::array<double, 256> _map = []() {
		std::array<double, 256> fs{};
		for (size_t v = 0; v < fs.size(); ++v) {
			fs[v] = static_cast<double>(v);
		}
		return fs;
	}();

	static constexpr bool hasRange = false;
	static double transform(uint16_t val) noexcept { return _map[val]; }
};

struct TLogitTransform {
	static constexpr bool hasRange = true;
	static double transform(uint16_t val) noexcept {
		// logit(coretools::Probability(genometools::PhredIntProbability(v)));
	constexpr std::array<double, 256> map = {
		5.8715919871401411001e+01,  1.3512152439604214749e+00,  5.3632602549409980064e-01,  4.7489433250102125808e-03,
		-4.1335816330680280606e-01, -7.7116213850872239455e-01, -1.0922828683340868317e+00, -1.3892580491035142476e+00,
		-1.6695115018271398100e+00, -1.9377746232188206577e+00, -2.1972245773362191201e+00, -2.4500782996574961281e+00,
		-2.6979279383599181763e+00, -2.9419423467915835069e+00, -3.1829942858483692980e+00, -3.4217440655229678370e+00,
		-3.6586964210559798083e+00, -3.8942402933127304721e+00, -4.1286772984382107410e+00, -4.3622425063144678603e+00,
		-4.5951198501345897895e+00, -4.8274536970397967650e+00, -5.0593576416241985427e+00, -5.2909212400500091888e+00,
		-5.5222152058020261833e+00, -5.7532954443123722754e+00, -5.9842061950776050949e+00, -6.2149824955657440029e+00,
		-6.4456521102104522569e+00, -6.6762370509949207076e+00, -6.9067547786485539163e+00, -7.1372191444315760123e+00,
		-7.3676411410677005165e+00, -7.5980294940648649415e+00, -7.8283911296266479596e+00, -8.0587315477558618682e+00,
		-8.2890551143869739548e+00, -8.5193652979228655653e+00, -8.7496648517883564011e+00, -8.9799559620469491250e+00,
		-9.2102403669758494686e+00, -9.4405194453277321287e+00, -9.6707942928185310905e+00, -9.9010657799494659059e+00,
		-1.0131334597547462906e+01, -1.0361601295249840859e+01, -1.0591866308397566598e+01, -1.0822129984234210909e+01,
		-1.1052392597604846713e+01, -1.1282654366173959559e+01, -1.1512915464920228104e+01, -1.1743176030986306912e+01,
		-1.1973436173944294936e+01, -1.2203695981037961715e+01, -1.2433955520971407438e+01, -1.2664214849237836802e+01,
		-1.2894474008681829247e+01, -1.3124733034786157049e+01, -1.3354991954761960216e+01, -1.3585250789575166408e+01,
		-1.3815509557963773446e+01, -1.4045768272965696966e+01, -1.4276026945574148641e+01, -1.4506285584729543103e+01,
		-1.4736544196937833462e+01, -1.4966802788286727832e+01, -1.5197061362376787841e+01, -1.5427319923518258094e+01,
		-1.5657578474161121918e+01, -1.5887837015602867652e+01, -1.6118095550958315698e+01, -1.6348354080855465043e+01,
		-1.6578612606429999232e+01, -1.6808871128792223004e+01, -1.7039129648228410474e+01, -1.7269388165885811048e+01,
		-1.7499646681440644613e+01, -1.7729905196085930186e+01, -1.7960163709795565978e+01, -1.8190422221900206523e+01,
		-1.8420680733952366381e+01, -1.8650939245339056782e+01, -1.8881197756210209349e+01, -1.9111456266893121381e+01,
		-1.9341714777052104068e+01, -1.9571973287340355796e+01, -1.9802231797041667249e+01, -2.0032490307037338795e+01,
		-2.0262748817053651607e+01, -2.0493007326224581988e+01, -2.0723265835946410363e+01, -2.0953524345482055224e+01,
		-2.1183782854882871050e+01, -2.1414041364397853329e+01, -2.1644299873629112341e+01, -2.1874558383180453092e+01,
		-2.2104816892296408781e+01, -2.2335075401827118924e+01, -2.2565333911474102280e+01, -2.2795592420351660223e+01,
		-2.3025850929840455450e+01, -2.3256109439190996824e+01, -2.3486367948444780041e+01, -2.3716626457842966857e+01,
		-2.3946884966981453857e+01, -2.4177143476459104221e+01, -2.4407401985516525400e+01, -2.4637660495000737626e+01,
		-2.4867919004610786970e+01, -2.5098177513459010157e+01, -2.5328436022924503135e+01, -2.5558694532256531318e+01,
		-2.5788953041495609853e+01, -2.6019211550882118900e+01, -2.6249470060011329764e+01, -2.6479728569481608247e+01,
		-2.6709987078533178106e+01, -2.6940245588012743383e+01, -2.7170504097619097905e+01, -2.7400762606464386550e+01,
		-2.7631021115927548948e+01, -2.7861279625257726167e+01, -2.8091538134495333878e+01, -2.8321796643880677635e+01,
		-2.8552055153008957689e+01, -2.8782313662478500760e+01, -2.9012572171529484422e+01, -2.9242830681008584293e+01,
		-2.9473089190614569333e+01, -2.9703347699459563103e+01, -2.9933606208922494574e+01, -3.0163864718252487052e+01,
		-3.0394123227489949102e+01, -3.0624381736875172066e+01, -3.0854640246003363302e+01, -3.1084898755472831766e+01,
		-3.1315157264523755032e+01, -3.1545415774002808718e+01, -3.1775674283608758230e+01, -3.2005932792453727131e+01,
		-3.2236191301916626628e+01, -3.2466449811246604895e+01, -3.2696708320484049182e+01, -3.2926966829869265041e+01,
		-3.3157225338997442066e+01, -3.3387483848466906977e+01, -3.3617742357517826690e+01, -3.3848000866996869718e+01,
		-3.4078259376602815678e+01, -3.4308517885447777473e+01, -3.4538776394910684076e+01, -3.4769034904240655237e+01,
		-3.4999293413478099524e+01, -3.5229551922863315383e+01, -3.5459810431991492408e+01, -3.5690068941460957319e+01,
		-3.5920327450511869927e+01, -3.6150585959990920060e+01, -3.6380844469596866020e+01, -3.6611102978441827815e+01,
		-3.6841361487904727312e+01, -3.7071619997234705579e+01, -3.7301878506472149866e+01, -3.7532137015857358620e+01,
		-3.7762395524985542750e+01, -3.7992654034455000556e+01, -3.8222912543505920269e+01, -3.8453171052984963296e+01,
		-3.8683429562590909256e+01, -3.8913688071435871052e+01, -3.9143946580898777654e+01, -3.9374205090228748816e+01,
		-3.9604463599466193102e+01, -3.9834722108851408962e+01, -4.0064980617979585986e+01, -4.0295239127449043792e+01,
		-4.0525497636499963505e+01, -4.0755756145979013638e+01, -4.0986014655584959598e+01, -4.1216273164429921394e+01,
		-4.1446531673892820891e+01, -4.1676790183222792052e+01, -4.1907048692460236339e+01, -4.2137307201845452198e+01,
		-4.2367565710973629223e+01, -4.2597824220443094134e+01, -4.2828082729494006742e+01, -4.3058341238973056875e+01,
		-4.3288599748579002835e+01, -4.3518858257423964631e+01, -4.3749116766886871233e+01, -4.3979375276216842394e+01,
		-4.4209633785454286681e+01, -4.4439892294839495435e+01, -4.4670150803967679565e+01, -4.4900409313437137371e+01,
		-4.5130667822488057084e+01, -4.5360926331967100111e+01, -4.5591184841573046072e+01, -4.5821443350418007867e+01,
		-4.6051701859880914469e+01, -4.6281960369210885631e+01, -4.6512218878448329917e+01, -4.6742477387833545777e+01,
		-4.6972735896961722801e+01, -4.7202994406431180607e+01, -4.7433252915482100320e+01, -4.7663511424961150453e+01,
		-4.7893769934567096414e+01, -4.8124028443412058209e+01, -4.8354286952874957706e+01, -4.8584545462204935973e+01,
		-4.8814803971442380259e+01, -4.9045062480827589013e+01, -4.9275320989955766038e+01, -4.9505579499425230949e+01,
		-4.9735838008476150662e+01, -4.9966096517955193690e+01, -5.0196355027561139650e+01, -5.0426613536406101446e+01,
		-5.0656872045869008048e+01, -5.0887130555198979209e+01, -5.1117389064436423496e+01, -5.1347647573821632250e+01,
		-5.1577906082949816380e+01, -5.1808164592419274186e+01, -5.2038423101470193899e+01, -5.2268681610949236926e+01,
		-5.2498940120555182887e+01, -5.2729198629400144682e+01, -5.2959457138863051284e+01, -5.3189715648193022446e+01,
		-5.3419974157430466732e+01, -5.3650232666815682592e+01, -5.3880491175943859616e+01, -5.4110749685413317422e+01,
		-5.4341008194464237135e+01, -5.4571266703943287268e+01, -5.4801525213549233229e+01, -5.5031783722394195024e+01,
		-5.5262042231857094521e+01, -5.5492300741187072788e+01, -5.5722559250424517074e+01, -5.5952817759809725828e+01,
		-5.6183076268937902853e+01, -5.6413334778407367764e+01, -5.6643593287458287477e+01, -5.6873851796937330505e+01,
		-5.7104110306543276465e+01, -5.7334368815388238261e+01, -5.7564627324851144863e+01, -5.7794885834181116024e+01,
		-5.8025144343418560311e+01, -5.8255402852803769065e+01, -5.8485661361931953195e+01, -5.8715919871401411001e+01};
		assert(val < 256);
		return map[val];
	}
};
} // namespace impl

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
template<size_t O, typename Transformer = impl::TLogitTransform> class TPolynomial : public TFunction {
	static_assert(O > 0);

private:
	std::array<double, O> _betas{1.};  // betas of the model
	std::array<double, O> _oldBetas{}; // use during estimation

	double _getAsDouble(uint16_t val) const noexcept {
		return Transformer::transform(val);
	}

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + O; }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + O; }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + O; }

public:
	static inline const std::string name = "polynomial";

	TPolynomial(uint16_t FirstParameterIndex, const std::vector<std::string> &betas) : TFunction(FirstParameterIndex) {
		_initializeValues(betas);
	}

	uint16_t numParameters() const noexcept override { return O; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return O; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) noexcept override {
		if constexpr (Transformer::hasRange) {
			return std::all_of(values.begin(), values.end(),
							   [](auto v) { return v <= genometools::PhredIntProbability::max().get(); });
		} else
			return true;
	};

	double adjustParametersPostEstimation() noexcept override { return 0.; }

	double getEtaTerm(uint16_t val) const noexcept override {
		const double v = _getAsDouble(val);
		if constexpr (O == 1) {
			return _betas.front() * v;
		} else if constexpr (O == 2) {
			return v * (_betas[0] + v * _betas[1]);
		} else if constexpr (O == 3) {
			return v * (_betas[0] + v * (_betas[1] + _betas[2] * v));
		} else {
			double vpi = v;
			double sum = _betas[0] * vpi;

			for (size_t i = 1; i < O; ++i) {
				vpi *= v;
				sum += _betas[i] * vpi;
			}
			return sum;
		}
	};

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override {
		const double v = _getAsDouble(val);
		if constexpr (O == 1) { return {firstParameterIndex(), v}; }
		else return {firstParameterIndex() + i, coretools::uPow(v, i + 1)};
	}

	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

	std::string typeString() const noexcept override {
		using coretools::str::toString;
		return name + '(' + toString(O) + ')';
	}
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
// A polynomial function
//--------------------------------------------------------------

class TProbit : public TFunction {
private:
	struct TProbitTmpStorage {
		double cumulDens_Phi;
		double normalDens_phi;
		double eta;
		double normalDens_q;
		double normalDens_Beta1;
		double normalDens_Beta1_q;
		double normalDens_Beta1_z;
		double normalDens_Beta1_q_z;
		double normalDens_Beta1_q2_z;

		TProbitTmpStorage(const std::array<double, 3> &betas, uint16_t q);
	};
	std::array<double, 3> _betas;    // betas of the model
	std::array<double, 3> _oldBetas; // use during estimation

	// tmp storage
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _expandTmpStorage(uint16_t MaxValue) const;

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _oldBetas.size(); }

public:
	static inline const std::string name = "probit";
	TProbit(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return 3; };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 3; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 6; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;
	T2ndDerivative get2ndDerivatives(uint16_t val, size_t i, size_t j) const noexcept override;

	bool checkOrInitValueRange(const std::vector<uint16_t> &) noexcept override { return true; };

	double getEtaTerm(uint16_t val) const noexcept override;
	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specific
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------
class TSpecific : public TFunction {
private:
	std::vector<double> _betas;    // betas of the model
	std::vector<double> _oldBetas; // use during estimation

	void _resize(uint16_t MaxValue);
	bool _checkValueRange(uint16_t val) const noexcept { return val < numParameters(); };
	void _adjustValueRanges(const std::vector<uint16_t> &values);

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _betas.size(); }

public:
	static inline const std::string name = "specific";

	TSpecific(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;
	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) override {
		if (numParameters() == 0) {
			_adjustValueRanges(values);
			return true;
		}

		for (uint16_t val : values) {
			if (!_checkValueRange(val)) return false;
		}
		return true;
	}

	double adjustParametersPostEstimation() noexcept override;

	double getEtaTerm(uint16_t val) const noexcept override { return _betas[val]; };
	virtual std::string typeString() const noexcept override { return name; }
};

//--------------------------------------------------------------
// TCovariateFunction_specificMap
// A term per discrete values as indicated with a map
//--------------------------------------------------------------
struct TIndexMapEntry {
	uint16_t index = 0;
	bool used      = false;
};

class TSpecificMap : public TFunction {
private:
	std::vector<double> _betas;            // betas of the model
	std::vector<double> _oldBetas;         // use during estimation
	std::vector<TIndexMapEntry> _indexMap; // maps value to parameter index

	void _resize(size_t NumParameters);
	void _initMapFromVector(const std::vector<uint16_t> &values);
	bool _checkValueRange(uint16_t val) const noexcept { return val < _indexMap.size() && _indexMap[val].used; };
	void _adjustValueRanges(const std::vector<uint16_t> &values);

protected:
	double *_begin() noexcept override { return _betas.data(); }
	double *_end() noexcept override { return _betas.data() + _betas.size(); }
	const double *_cbegin() const noexcept override { return _betas.data(); }
	const double *_cend() const noexcept override { return _betas.data() + _betas.size(); }
	double *_obegin() noexcept override { return _oldBetas.data(); }
	double *_oend() noexcept override { return _oldBetas.data() + _betas.size(); }

public:
	static inline const std::string name = "map";
	TSpecificMap(uint16_t FirstParameterIndex, const std::vector<std::string> &betas);

	uint16_t numParameters() const noexcept override { return _betas.size(); };
	uint16_t numNonZeroFirstDerivatives() const noexcept override { return 1; };
	uint16_t numNonZeroSecondDerivatives() const noexcept override { return 0; };

	T1stDerivative get1stDerivatives(uint16_t val, size_t i) const noexcept override;

	T2ndDerivative get2ndDerivatives(uint16_t, size_t, size_t ) const noexcept override {
		return {firstParameterIndex(), firstParameterIndex(), 0.};
	};

	bool checkOrInitValueRange(const std::vector<uint16_t> &values) override {
		if (numParameters() == 0) {
			_adjustValueRanges(values);
			return true;
		}
		for (uint16_t val : values) {
			if (!_checkValueRange(val)) return false;
		}
		return true;
	}

	double adjustParametersPostEstimation() noexcept override;

	double getEtaTerm(uint16_t val) const noexcept override { return _betas[_indexMap[val].index]; };
	virtual std::string typeString() const noexcept override { return name; }
};

} // namespace SequencingError
} // namespace GenotypeLikelihoods

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
