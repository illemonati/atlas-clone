/*
 * Types.h
 *
 *  Created on: May 7, 2021
 *      Author: phaentu
 */

#ifndef BAM_TYPES_H_
#define BAM_TYPES_H_

#include <cinttypes>
#include <ostream>
#include <math.h>
#include <cmath>
#include "strongTypes.h"

//------------------------------------------------
// Base
//------------------------------------------------

enum BaseEnum : uint8_t {A=0, C, G, T, N};

class Base : public StrongTypes::StrongType<BaseEnum, Base> {
private:
	static constexpr char _toChar[5] = {'A', 'C', 'G', 'T', 'N'};
	static constexpr BaseEnum _flipBase[5] = {T, G, C, A, N};
public:
	//constructors
	explicit constexpr Base() : StrongType(N) {};
	explicit constexpr Base(const BaseEnum & base) : StrongType(base) {};
	explicit constexpr Base(const char & base){
		if(base == 'A' || base == 'a') _value = A;
		else if(base == 'C' || base == 'c') _value = C;
		else if(base == 'G' || base == 'g') _value = G;
		else if(base == 'T' || base == 't') _value = T;
		else _value = N;
	};

	//convert
	explicit constexpr operator char() const {
		return _toChar[ static_cast<uint8_t>(_value) ];
	};

	//manipulate
	void flip(){
		_value = _flipBase[ static_cast<uint8_t>(_value) ];
	};
	Base constexpr flipped() const{
		return Base(_flipBase[ static_cast<uint8_t>(_value) ]);
	};
};

std::ostream& operator<<(std::ostream& os, const Base & base);

//------------------------------------------------
// Genotype
//------------------------------------------------

enum GenotypeEnum : uint8_t {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT, NN};

class Genotype : public StrongTypes::StrongType<GenotypeEnum, Genotype> {
private:

	static std::string _toString[NN+1];

	constexpr GenotypeEnum _toGeno(const Base & first, const Base & second){
		//assume first <= second, e.g. AC but not CA
		if(first == N || second == N){
			return NN;
		} else if(first == A){
			return static_cast<GenotypeEnum>(     static_cast<uint8_t>((BaseEnum) second) );
		} else if(first == C){
			return static_cast<GenotypeEnum>( 3 + static_cast<uint8_t>((BaseEnum) second) );
		} else if(first == G){
			return static_cast<GenotypeEnum>( 5 + static_cast<uint8_t>((BaseEnum) second) );
		} else { //first == T
			return static_cast<GenotypeEnum>( 6 + static_cast<uint8_t>((BaseEnum) second) );
		}
	};

public:
	//constructors
	explicit constexpr Genotype() : StrongType(NN) {};

	explicit constexpr Genotype(const GenotypeEnum & genotype) : StrongType(genotype) {};

	explicit constexpr Genotype(const Base & first, const Base & second){
		if(static_cast<uint8_t>( (BaseEnum) first) > static_cast<uint8_t>( (BaseEnum) second)){
			_value = _toGeno(second, first);
		} else {
			_value = _toGeno(first, second);
		}
	};

	explicit constexpr Genotype(const BaseEnum & first, const BaseEnum & second) : Genotype(Base(first), Base(second)) {};

	//convert
	explicit operator std::string() const {
		return _toString[static_cast<uint8_t>( (GenotypeEnum) _value)];
	};
};

std::ostream& operator<<(std::ostream& os, const Genotype & genotype);

//------------------------------------------------
// Context
//------------------------------------------------

enum BaseContextEnum : uint8_t {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT, cNN}; //N means unknown base or "nothing", i.e. end of read

class BaseContext : public StrongTypes::StrongType<BaseContextEnum, BaseContext> {
private:
	static const std::string _toString[cNN+1];

public:
	//constructors
	explicit constexpr BaseContext() : StrongType(cNN) {};

	explicit constexpr BaseContext(const BaseContextEnum & context) : StrongType(context) {};

	explicit constexpr BaseContext(const Base & first, const Base & second){
		if(second == N){
			_value = cNN;
		} else {
			_value = static_cast<BaseContextEnum>( static_cast<uint8_t>(N) * static_cast<uint8_t>( (BaseEnum) first) + static_cast<uint8_t>( (BaseEnum) second));
		}
	};

	explicit constexpr BaseContext(const BaseEnum & first, const BaseEnum & second) : BaseContext(Base(first), Base(second)) {};

	//convert
	explicit operator std::string() const {
		return _toString[static_cast<uint8_t>( (BaseContextEnum) _value)];
	};
};

std::ostream& operator<<(std::ostream& os, const BaseContext & context);

//-------------------------------------
// ErrorRate
// An error rate within [0,1]
//-------------------------------------
class LogErrorRate;
class PhredErrorRate;
class PhredIntErrorRate;
class Quality;

class ErrorRate : public StrongTypes::StrongType<double, ErrorRate, StrongTypes::Printable>{
private:
	static const double _phredIntToError[256];

public:
	//constructors
	explicit constexpr ErrorRate() : StrongType() {};

	explicit constexpr ErrorRate(const double & error){
		if(error < 0.0 || error > 1.0){
			throw std::runtime_error("explicit constexpr ErrorRate(const double & error): error is outside [0,1]!");
		}
		_value = error;
	};

	explicit ErrorRate(const LogErrorRate & error);
	explicit ErrorRate(const PhredErrorRate & error);
	explicit ErrorRate(const PhredIntErrorRate & error);
	explicit ErrorRate(const Quality & quality);
};

//-------------------------------------
// LogErrorRate
// The log of an error rate
//-------------------------------------
class LogErrorRate : public StrongTypes::StrongType<double, LogErrorRate, StrongTypes::Printable>{
private:
	static const double _phredIntToLogError[256];

public:
	//constructors
	explicit constexpr LogErrorRate() : StrongType() {};

	explicit constexpr LogErrorRate(const double & logError){
		if(logError > 0.0){
			throw std::runtime_error("explicit constexpr ErrorRate(const double & error): error is > 0!");
		}
		_value = logError;
	};

	explicit constexpr LogErrorRate(const ErrorRate & error){
		_value = log( (double) error );
	};

	explicit LogErrorRate(const PhredErrorRate & error);
	explicit LogErrorRate(const PhredIntErrorRate & error);
	explicit LogErrorRate(const Quality & quality);
};

//------------------------------------------------
// PhredErrorRate
// phreded error = -10 * log_10(error)
//------------------------------------------------
class PhredErrorRate : public StrongTypes::StrongType<double, PhredErrorRate, StrongTypes::Printable> {
private:

public:
	//constructors
	explicit constexpr PhredErrorRate() : StrongType() {};
	explicit constexpr PhredErrorRate(const double & phredError){
		if(phredError < 0){
			throw std::runtime_error("explicit constexpr PhredErrorRate(const double & phredError): phredError < 0!");
		}
		_value = phredError;
	};
	explicit constexpr PhredErrorRate(const ErrorRate & error){
		_value = -10.0 * log10( (double) error);
	};

	explicit constexpr PhredErrorRate(const LogErrorRate & logError){
		_value = -4.342944819032518 * (double) logError; // -10 * log_10(e) * log(error)
	};

	explicit PhredErrorRate(const PhredIntErrorRate & error);
	explicit PhredErrorRate(const Quality & quality);
};

//------------------------------------------------
// PhredIntError
// phreded error stored as uint8_t
// only valid within [0,255] and truncated outside
//------------------------------------------------
class PhredIntErrorRate : public StrongTypes::StrongType<uint8_t, PhredIntErrorRate, StrongTypes::Printable> {
private:

public:
	//constructors
	explicit constexpr PhredIntErrorRate() : StrongType() {};
	explicit constexpr PhredIntErrorRate(const uint8_t & phredInt) : StrongType(phredInt) {};
	explicit PhredIntErrorRate(const ErrorRate & error);
	explicit PhredIntErrorRate(const LogErrorRate & logError);

	explicit constexpr PhredIntErrorRate(const PhredErrorRate & error){
		if( (double) error > 255.0){
			_value = 255;
		} else {
			_value = (uint8_t) std::round( (double) error);
		}
	};

	explicit PhredIntErrorRate(const Quality & quality);
};

//------------------------------------------------
// Quality
// Sequencing or mapping quality = PhredIntError + 33
// only valid within [33,126] (matching chars ! until ~) and truncated outside
//------------------------------------------------
class Quality : public StrongTypes::StrongType<char, Quality, StrongTypes::Printable> {
private:

public:
	//constructors
	explicit constexpr Quality() : StrongType() {};
	explicit Quality(const ErrorRate & error);
	explicit Quality(const LogErrorRate & error);
	explicit Quality(const PhredErrorRate & error);

	explicit constexpr Quality(const PhredIntErrorRate & error){
		if( (uint8_t) error > 125){
			_value = 126;
		} else {
			_value = (uint8_t) error + 33;
		}
	};

	explicit constexpr Quality(const char & quality){
		if(quality < 33){
			_value = 33;
		} else if(quality > 125){
			_value = 126;
		} else {
			_value = quality;
		}
	};
};




#endif /* BAM_TYPES_H_ */
