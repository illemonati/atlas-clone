/*
 * TRecalibrationEMModule.h
 *
 *  Created on: Jan 20, 2020
 *      Author: wegmannd
 */

#ifndef RECALIBRATION_TRECALIBRATIONEMMODULE_H_
#define RECALIBRATION_TRECALIBRATIONEMMODULE_H_

//--------------------------------------------------------------
// TRecalibrationEMModule
// Base class for recal modules
//--------------------------------------------------------------

class TRecalibrationEMModule{
protected:
	int _numParameters;
	double* _betas; //betas of the model
	double* _oldBetas; //use during estimation
	bool _initialized;

	void _freeBetas();
	void _initializeBetas();

	friend class TRecalibrationEMModelNEW;

public:
	TRecalibrationEMModule(){
		_initialized = false;
		_numParameters = 0;
		_betas = nullptr;
		_oldBetas = nullptr;
	};

	virtual ~TRecalibrationEMModule(){
		_freeBetas();
	};

	virtual double getEtaTerm(const uint8_t val){ return 0.0; }
	virtual double getEtaTerm(const uint16_t val){ return 0.0; }
};


//--------------------------------------------------------------
// TRecalibrationEMModule_intercept
// An intercept term
//--------------------------------------------------------------
class TRecalibrationEMModule_intercept:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_intercept(){
		_numParameters = 1;
		_initializeBetas();
	};

	double getEtaTerm(){
		return _betas[0];
	};
	double getEtaTerm(const uint8_t val){
		return _betas[0];
	};
	double getEtaTerm(const uint16_t val){
		return _betas[0];
	};
};


//--------------------------------------------------------------
// TRecalibrationEMModule_linear
// A linear function
//--------------------------------------------------------------
class TRecalibrationEMModule_linear:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_linear(){
		_numParameters = 1;
		_initializeBetas();
	};

	double getEtaTerm(const uint8_t val){
		return _betas[0] * val;
	};
	double getEtaTerm(const uint16_t val){
		return _betas[0] * val;
	};
};


//--------------------------------------------------------------
// TRecalibrationEMModule_quadratic
// A quadratic function
//--------------------------------------------------------------
class TRecalibrationEMModule_quadratic:public TRecalibrationEMModule{
protected:

public:
	TRecalibrationEMModule_quadratic(){
		_numParameters = 2;
		_initializeBetas();
	};

	double getEtaTerm(const uint8_t val){
		return val * (_betas[0] + _betas[1] * val);
	};
	double getEtaTerm(const uint16_t val){
		return val * (_betas[0] + _betas[1] * val);
	};
};

//--------------------------------------------------------------
// TRecalibrationEMModule_specific
// A term per discrete value from o to maxValue
//--------------------------------------------------------------
class TRecalibrationEMModule_specific:public TRecalibrationEMModule{
protected:
	int _maxValue;

public:
	TRecalibrationEMModule_specific(int maxValue){
		_maxValue = maxValue;
		_numParameters = _maxValue + 1;
		_initializeBetas();
	};

	double getEtaTerm(const uint8_t val){
		return _betas[val];
	};
	double getEtaTerm(const uint16_t val){
		return _betas[val];
	};
};

#endif /* RECALIBRATION_TRECALIBRATIONEMMODULE_H_ */
