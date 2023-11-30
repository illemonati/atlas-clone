/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_
#define GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_

#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "coretools/Main/TLog.h"
#include "coretools/Math/TProbabilityDistributions.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/toString.h"
#include "coretools/Types/probability.h"
#include "genometools/PhredProbabilityTypes.h"

#include "SequencingError/RecalEstimatorTools.h"
#include "SequencingError/TCovariate.h"
#include "TDerivatives.h"
#include "TReadGroupInfo.h"

namespace GenotypeLikelihoods::SequencingError {

//--------------------------------------------------------------
// TCovariateFunction
// Base class for recal covariate functions
//--------------------------------------------------------------
class TFunction {
protected:
	size_t _firstParameterIndex;

public:
	TFunction(size_t FirstParameterIndex) : _firstParameterIndex(FirstParameterIndex) {}
	virtual ~TFunction() = default;

	// non-virtuals
	constexpr size_t firstParameterIndex() const noexcept { return _firstParameterIndex; }

	// virtuals
	virtual double *begin() noexcept              = 0;
	virtual double *end() noexcept                = 0;
	virtual const double *begin() const noexcept  = 0;
	virtual const double *end() const noexcept    = 0;
	virtual size_t numParameters() const noexcept = 0;

	// check value range: to ensure that data can be recalibrated
	virtual void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex) = 0;

	// estimation
	virtual double getEta(const BAM::TSequencedBase &base) const noexcept   = 0;
	virtual double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
						  std::vector<T2ndDerivative> &der2) const noexcept = 0;
	virtual double adjust() noexcept                                        = 0;
	virtual std::string typeString() const noexcept                         = 0;
	virtual void addInfo(BAM::RGInfo::TInfo &info) const                    = 0;
	virtual void log() const; 
};

class TNoFunction final : public TFunction{
	double _nothing = 0;
public:
	TNoFunction() : TFunction(0) {}

	// virtuals
	double *begin() noexcept override {return &_nothing;}
	double *end() noexcept override {return &_nothing;}
	const double *begin() const noexcept override {return &_nothing;}
	const double *end() const noexcept override {return &_nothing;}
	virtual size_t numParameters() const noexcept override {return 0;}

	// check value range: to ensure that data can be recalibrated
	void init(const RecalEstimatorTools::TRecalDataTable &, size_t) override {}

	// estimation
	double getEta(const BAM::TSequencedBase &) const noexcept override {return 0.;}
	double getEta(const BAM::TSequencedBase &, std::vector<T1stDerivative> &,
						  std::vector<T2ndDerivative> &) const noexcept override {return 0.;}
	double adjust() noexcept override {return 0.;}
	std::string typeString() const noexcept override {return "";}
	void log() const noexcept  override {}; 
	void addInfo(BAM::RGInfo::TInfo &) const override {};
};


//--------------------------------------------------------------
// TCovariateFunction_intercept
// An intercept term
//--------------------------------------------------------------
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

	void addInfo(BAM::RGInfo::TInfo& info) const {
		 info[name] = _beta;
	}

	void log() const {
		coretools::instances::logfile().list(typeString(), ": ", _beta);
	}
};

namespace impl {
struct TNoTransform {
	static constexpr size_t N = uint16_t(-1);
	static double transform(uint16_t val) noexcept { return static_cast<double>(val); }
};

struct TLogitTransform {
	static constexpr size_t N = genometools::PhredIntProbability::max().get();
	static double transform(uint16_t val) noexcept {
		// logit(coretools::Probability(genometools::PhredIntProbability(v)));
		constexpr std::array<double, 256> map = {
			5.8715919871401411001e+01,  1.3512152439604214749e+00,  5.3632602549409980064e-01,
			4.7489433250102125808e-03,  -4.1335816330680280606e-01, -7.7116213850872239455e-01,
			-1.0922828683340868317e+00, -1.3892580491035142476e+00, -1.6695115018271398100e+00,
			-1.9377746232188206577e+00, -2.1972245773362191201e+00, -2.4500782996574961281e+00,
			-2.6979279383599181763e+00, -2.9419423467915835069e+00, -3.1829942858483692980e+00,
			-3.4217440655229678370e+00, -3.6586964210559798083e+00, -3.8942402933127304721e+00,
			-4.1286772984382107410e+00, -4.3622425063144678603e+00, -4.5951198501345897895e+00,
			-4.8274536970397967650e+00, -5.0593576416241985427e+00, -5.2909212400500091888e+00,
			-5.5222152058020261833e+00, -5.7532954443123722754e+00, -5.9842061950776050949e+00,
			-6.2149824955657440029e+00, -6.4456521102104522569e+00, -6.6762370509949207076e+00,
			-6.9067547786485539163e+00, -7.1372191444315760123e+00, -7.3676411410677005165e+00,
			-7.5980294940648649415e+00, -7.8283911296266479596e+00, -8.0587315477558618682e+00,
			-8.2890551143869739548e+00, -8.5193652979228655653e+00, -8.7496648517883564011e+00,
			-8.9799559620469491250e+00, -9.2102403669758494686e+00, -9.4405194453277321287e+00,
			-9.6707942928185310905e+00, -9.9010657799494659059e+00, -1.0131334597547462906e+01,
			-1.0361601295249840859e+01, -1.0591866308397566598e+01, -1.0822129984234210909e+01,
			-1.1052392597604846713e+01, -1.1282654366173959559e+01, -1.1512915464920228104e+01,
			-1.1743176030986306912e+01, -1.1973436173944294936e+01, -1.2203695981037961715e+01,
			-1.2433955520971407438e+01, -1.2664214849237836802e+01, -1.2894474008681829247e+01,
			-1.3124733034786157049e+01, -1.3354991954761960216e+01, -1.3585250789575166408e+01,
			-1.3815509557963773446e+01, -1.4045768272965696966e+01, -1.4276026945574148641e+01,
			-1.4506285584729543103e+01, -1.4736544196937833462e+01, -1.4966802788286727832e+01,
			-1.5197061362376787841e+01, -1.5427319923518258094e+01, -1.5657578474161121918e+01,
			-1.5887837015602867652e+01, -1.6118095550958315698e+01, -1.6348354080855465043e+01,
			-1.6578612606429999232e+01, -1.6808871128792223004e+01, -1.7039129648228410474e+01,
			-1.7269388165885811048e+01, -1.7499646681440644613e+01, -1.7729905196085930186e+01,
			-1.7960163709795565978e+01, -1.8190422221900206523e+01, -1.8420680733952366381e+01,
			-1.8650939245339056782e+01, -1.8881197756210209349e+01, -1.9111456266893121381e+01,
			-1.9341714777052104068e+01, -1.9571973287340355796e+01, -1.9802231797041667249e+01,
			-2.0032490307037338795e+01, -2.0262748817053651607e+01, -2.0493007326224581988e+01,
			-2.0723265835946410363e+01, -2.0953524345482055224e+01, -2.1183782854882871050e+01,
			-2.1414041364397853329e+01, -2.1644299873629112341e+01, -2.1874558383180453092e+01,
			-2.2104816892296408781e+01, -2.2335075401827118924e+01, -2.2565333911474102280e+01,
			-2.2795592420351660223e+01, -2.3025850929840455450e+01, -2.3256109439190996824e+01,
			-2.3486367948444780041e+01, -2.3716626457842966857e+01, -2.3946884966981453857e+01,
			-2.4177143476459104221e+01, -2.4407401985516525400e+01, -2.4637660495000737626e+01,
			-2.4867919004610786970e+01, -2.5098177513459010157e+01, -2.5328436022924503135e+01,
			-2.5558694532256531318e+01, -2.5788953041495609853e+01, -2.6019211550882118900e+01,
			-2.6249470060011329764e+01, -2.6479728569481608247e+01, -2.6709987078533178106e+01,
			-2.6940245588012743383e+01, -2.7170504097619097905e+01, -2.7400762606464386550e+01,
			-2.7631021115927548948e+01, -2.7861279625257726167e+01, -2.8091538134495333878e+01,
			-2.8321796643880677635e+01, -2.8552055153008957689e+01, -2.8782313662478500760e+01,
			-2.9012572171529484422e+01, -2.9242830681008584293e+01, -2.9473089190614569333e+01,
			-2.9703347699459563103e+01, -2.9933606208922494574e+01, -3.0163864718252487052e+01,
			-3.0394123227489949102e+01, -3.0624381736875172066e+01, -3.0854640246003363302e+01,
			-3.1084898755472831766e+01, -3.1315157264523755032e+01, -3.1545415774002808718e+01,
			-3.1775674283608758230e+01, -3.2005932792453727131e+01, -3.2236191301916626628e+01,
			-3.2466449811246604895e+01, -3.2696708320484049182e+01, -3.2926966829869265041e+01,
			-3.3157225338997442066e+01, -3.3387483848466906977e+01, -3.3617742357517826690e+01,
			-3.3848000866996869718e+01, -3.4078259376602815678e+01, -3.4308517885447777473e+01,
			-3.4538776394910684076e+01, -3.4769034904240655237e+01, -3.4999293413478099524e+01,
			-3.5229551922863315383e+01, -3.5459810431991492408e+01, -3.5690068941460957319e+01,
			-3.5920327450511869927e+01, -3.6150585959990920060e+01, -3.6380844469596866020e+01,
			-3.6611102978441827815e+01, -3.6841361487904727312e+01, -3.7071619997234705579e+01,
			-3.7301878506472149866e+01, -3.7532137015857358620e+01, -3.7762395524985542750e+01,
			-3.7992654034455000556e+01, -3.8222912543505920269e+01, -3.8453171052984963296e+01,
			-3.8683429562590909256e+01, -3.8913688071435871052e+01, -3.9143946580898777654e+01,
			-3.9374205090228748816e+01, -3.9604463599466193102e+01, -3.9834722108851408962e+01,
			-4.0064980617979585986e+01, -4.0295239127449043792e+01, -4.0525497636499963505e+01,
			-4.0755756145979013638e+01, -4.0986014655584959598e+01, -4.1216273164429921394e+01,
			-4.1446531673892820891e+01, -4.1676790183222792052e+01, -4.1907048692460236339e+01,
			-4.2137307201845452198e+01, -4.2367565710973629223e+01, -4.2597824220443094134e+01,
			-4.2828082729494006742e+01, -4.3058341238973056875e+01, -4.3288599748579002835e+01,
			-4.3518858257423964631e+01, -4.3749116766886871233e+01, -4.3979375276216842394e+01,
			-4.4209633785454286681e+01, -4.4439892294839495435e+01, -4.4670150803967679565e+01,
			-4.4900409313437137371e+01, -4.5130667822488057084e+01, -4.5360926331967100111e+01,
			-4.5591184841573046072e+01, -4.5821443350418007867e+01, -4.6051701859880914469e+01,
			-4.6281960369210885631e+01, -4.6512218878448329917e+01, -4.6742477387833545777e+01,
			-4.6972735896961722801e+01, -4.7202994406431180607e+01, -4.7433252915482100320e+01,
			-4.7663511424961150453e+01, -4.7893769934567096414e+01, -4.8124028443412058209e+01,
			-4.8354286952874957706e+01, -4.8584545462204935973e+01, -4.8814803971442380259e+01,
			-4.9045062480827589013e+01, -4.9275320989955766038e+01, -4.9505579499425230949e+01,
			-4.9735838008476150662e+01, -4.9966096517955193690e+01, -5.0196355027561139650e+01,
			-5.0426613536406101446e+01, -5.0656872045869008048e+01, -5.0887130555198979209e+01,
			-5.1117389064436423496e+01, -5.1347647573821632250e+01, -5.1577906082949816380e+01,
			-5.1808164592419274186e+01, -5.2038423101470193899e+01, -5.2268681610949236926e+01,
			-5.2498940120555182887e+01, -5.2729198629400144682e+01, -5.2959457138863051284e+01,
			-5.3189715648193022446e+01, -5.3419974157430466732e+01, -5.3650232666815682592e+01,
			-5.3880491175943859616e+01, -5.4110749685413317422e+01, -5.4341008194464237135e+01,
			-5.4571266703943287268e+01, -5.4801525213549233229e+01, -5.5031783722394195024e+01,
			-5.5262042231857094521e+01, -5.5492300741187072788e+01, -5.5722559250424517074e+01,
			-5.5952817759809725828e+01, -5.6183076268937902853e+01, -5.6413334778407367764e+01,
			-5.6643593287458287477e+01, -5.6873851796937330505e+01, -5.7104110306543276465e+01,
			-5.7334368815388238261e+01, -5.7564627324851144863e+01, -5.7794885834181116024e+01,
			-5.8025144343418560311e+01, -5.8255402852803769065e+01, -5.8485661361931953195e+01,
			-5.8715919871401411001e+01};
		assert(val < 256);
		return map[val];
	}
};
} // namespace impl

//--------------------------------------------------------------
// TCovariateFunction_polynomial
// A polynomial function
//--------------------------------------------------------------
template<size_t O, typename Covariate> class TPolynomial final : public TFunction {
	static_assert(O > 0);

private:
	using Transformer =
		std::conditional_t<std::is_same_v<Covariate, TCovariate_quality>, impl::TLogitTransform, impl::TNoTransform>;
	std::array<double, O> _betas{};  // betas of the model

public:
	static constexpr std::string_view name = "polynomial";

	TPolynomial(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {
		if constexpr (std::is_same_v<Covariate, TCovariate_quality>) {
			_betas.front() = 1.;
		}
	}

	size_t numParameters() const noexcept override { return O; }

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + O; }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + O; }


	void init(const RecalEstimatorTools::TRecalDataTable &, size_t FirstParameterIndex) noexcept override {
		_firstParameterIndex = FirstParameterIndex;
	}

	double adjust() noexcept override { return 0.; }

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const double v = Transformer::transform(Covariate::extract(base));
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
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
						  std::vector<T2ndDerivative> &) const noexcept override {
		const double v = Transformer::transform(Covariate::extract(base));
		if constexpr (O == 1) {
			der1.emplace_back(firstParameterIndex(), v);
			return _betas.front() * v;
		} else if constexpr (O == 2) {
			der1.emplace_back(firstParameterIndex(), v);
			der1.emplace_back(firstParameterIndex() + 1, v*v);
			return v * (_betas[0] + v * _betas[1]);
		} else if constexpr (O == 3) {
			der1.emplace_back(firstParameterIndex(), v);
			der1.emplace_back(firstParameterIndex() + 1, v*v);
			der1.emplace_back(firstParameterIndex() + 2, v*v*v);
			return v * (_betas[0] + v * (_betas[1] + _betas[2] * v));
		} else {
			der1.emplace_back(firstParameterIndex(), v);
			double vpi = v;
			double sum = _betas[0] * vpi;

			for (size_t i = 1; i < O; ++i) {
				vpi *= v;
				sum += _betas[i] * vpi;
				der1.emplace_back(firstParameterIndex() + i, vpi);
			}
			return sum;
		}
	}

	std::string typeString() const noexcept override {
		using coretools::str::toString;
		return toString(Covariate::name, ":", name, O);
	}

	void addInfo(BAM::RGInfo::TInfo& info) const override {
		info[Covariate::name] = {{name, _betas}};
	}

	void log() const override {
		using coretools::instances::logfile;
		using coretools::str::toString;
		const auto var = TCovariate_quality::name.substr(0, 1);
		auto out = toString(_betas[0], "*", var);
		for (size_t i = 1; i < _betas.size(); ++i) {
			if (_betas[i] < 0) out += toString(" - ", -_betas[i], "*", var, "^", i+1);
			else out += toString(" + ", _betas[i], "*", var, "^", i+1);
		}
		logfile().list(typeString(), ": ", out);
	}
};

//--------------------------------------------------------------
// TRecalibrationEMCovariateFunction_probit
// A polynomial function
//--------------------------------------------------------------

template<typename Covariate> class TProbit final : public TFunction {
private:
	struct TProbitTmpStorage {
		double phi;
		double phiCumul;

		TProbitTmpStorage(const std::array<double, 3> &betas, uint16_t q) {
			using namespace coretools::probdist;
			const double z = betas[1] + betas[2] * q;
			phiCumul       = coretools::probdist::TNormalDistr::cumulativeDensity(z, 0, 1);
			phi            = coretools::probdist::TNormalDistr::density(z, 0, 1);
		}
	};
	std::array<double, 3> _betas{1., 0, 1.};    // betas of the model

	// tmp storage
	mutable std::vector<TProbitTmpStorage> _tmpStorage;

	void _expandTmpStorage(size_t MaxValue) const {
		for (size_t q = _tmpStorage.size(); q <= MaxValue; ++q) { _tmpStorage.emplace_back(_betas, q); }
	}

public:
	static constexpr std::string_view name = "probit";
	TProbit(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {}

	size_t numParameters() const noexcept override { return 3; }

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + _betas.size(); }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + _betas.size(); }

	void init(const RecalEstimatorTools::TRecalDataTable &, size_t FirstParameterIndex) noexcept override {
		_firstParameterIndex = FirstParameterIndex;
	}

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const auto q = Covariate::extract(base);
		if (q >= _tmpStorage.size()) { _expandTmpStorage(q); }
		return _tmpStorage[q].phiCumul*_betas.front();
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &der2) const noexcept override {
		const auto q = Covariate::extract(base);
		if (q >= _tmpStorage.size()) { _expandTmpStorage(q); }

		const auto b1 = firstParameterIndex();
		const auto b2 = b1 + 1;
		const auto b3 = b1 + 2;

		const auto z        = _betas[1] + _betas[2] * q;
		const auto phi      = _tmpStorage[q].phi;
		const auto phiCumul = _tmpStorage[q].phiCumul;
		const double bPhi   = phi * _betas.front();
		const double bPhiQ  = bPhi * q;
		const double bPhiQZ = bPhiQ * z;

		der1.emplace_back(b1, phiCumul);
		der1.emplace_back(b2, bPhi);
		der1.emplace_back(b3, bPhiQ);

		//der2.emplace_back(b1, b1, 0.); this is zero
		der2.emplace_back(b1, b2, phi);
		der2.emplace_back(b1, b3, phi*q);

		//der2.emplace_back(b2, b1, phi); only upper half
		der2.emplace_back(b2, b2, -bPhi*z);
		der2.emplace_back(b2, b3, -bPhiQZ*q);

		//der2.emplace_back(b3, b1, phi*p); only upper half
		//der2.emplace_back(b3, b2, -bPhi*z); only upper half
		der2.emplace_back(b3, b3, -bPhiQZ);

		return phiCumul *_betas.front();
	}

	double adjust() noexcept override { return 0.; }
	std::string typeString() const noexcept override { return std::string(Covariate::name).append(1, ':').append(name); }

	void addInfo(BAM::RGInfo::TInfo& info) const override {
		info[Covariate::name] = {{name, _betas}};
	}
};

//--------------------------------------------------------------
// A term per discrete value from 0 to maxValue
//--------------------------------------------------------------

template<typename Covariate> class TEmpiric final : public TFunction {
private:
	std::vector<double> _betas;    // betas of the model

public:
	static constexpr std::string_view name = "empiric";

	TEmpiric(size_t FirstParameterIndex) : TFunction(FirstParameterIndex) {}

	size_t numParameters() const noexcept override { return _betas.size(); }

	double *begin() noexcept override { return _betas.data(); }
	double *end() noexcept override { return _betas.data() + _betas.size(); }
	const double *begin() const noexcept override { return _betas.data(); }
	const double *end() const noexcept override { return _betas.data() + _betas.size(); }

	void push_back(double val) noexcept {_betas.push_back(val);}

	void init(const RecalEstimatorTools::TRecalDataTable &dataTable, size_t FirstParameterIndex) override {
		_firstParameterIndex = FirstParameterIndex;
		_betas.assign(Covariate::N(dataTable[Covariate::index]), NAN);
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (Covariate::isUsed(dataTable[Covariate::index], i)) {
				if constexpr (std::is_same_v<Covariate, TCovariate_quality>) {
					const coretools::Probability p = coretools::Probability(genometools::PhredIntProbability(i));
					_betas[i] = coretools::logit(p);
				} else {
					_betas[i] = 0.;
				}
			}
		}
	}

	double adjust() noexcept override {
		double mean = 0.;
		size_t N    = 0;
		for (auto bi : _betas) {
			if (!std::isnan(bi)) {
				++N;
				mean += bi;
			}
		}
		if (N > 1) mean /= N;

		for (auto &bi : _betas) {
			if (!std::isnan(bi)) { bi -= mean; }
		}
		return mean;
	}

	double getEta(const BAM::TSequencedBase &base) const noexcept override {
		const auto val = Covariate::extract(base);
		if (val < _betas.size()) return _betas[val];

		return _betas.back();
	}

	double getEta(const BAM::TSequencedBase &base, std::vector<T1stDerivative> &der1,
				  std::vector<T2ndDerivative> &) const noexcept override {
		const auto val = Covariate::extract(base);
		assert(val < _betas.size());

		const size_t der_index = firstParameterIndex() + static_cast<size_t>(val);
		der1.emplace_back(der_index, 1.0);
		return _betas[val];
	}

	std::string typeString() const noexcept override {
		return std::string(Covariate::name).append(1, ':').append(name);
	}

	void addInfo(BAM::RGInfo::TInfo &info) const override {
		BAM::RGInfo::TInfo ar = nlohmann::json::array();
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (!std::isnan(_betas[i])) ar += {i, _betas[i]};
		}
		info[Covariate::name] = {{name, ar}};
	}

	void log() const override {
		using coretools::str::toString;
		using coretools::instances::logfile;
		constexpr size_t Nmax = 3;

		std::vector<size_t> iis;
		iis.reserve(_betas.size());
		for (size_t i = 0; i < _betas.size(); ++i) {
			if (!std::isnan(_betas[i])) iis.push_back(i);
		}
		if (iis.empty()) {
			logfile().list(typeString(), ": []");
			return;
		}

		std::string ret = "[";
		if (iis.size() <= 2 * Nmax) {
			for (auto i : iis) ret.append(toString(i, ": ", _betas[i], ", "));
		} else {
			for (size_t j = 0; j < Nmax; ++j)
				ret.append(toString(iis[j], ": ", _betas[iis[j]], ", "));
			ret.append("..., ");
			const auto jStart = iis.size() - 1 - Nmax;
			for (size_t j = 0; j < Nmax; ++j)
				ret.append(toString(iis[jStart + j], ": ", _betas[iis[jStart + j]], ", "));
		}
		ret.pop_back();
		ret.back() = ']';
		logfile().list(typeString(), ": ", ret);
	}
};

} // namespace GenotypeLikelihoods::SequencingError

#endif /* GENOTYPELIKELIHOODS_TSEQUENCINGERRORCOVARIATEFUNCTION_H_ */
