/*
 * Types.h
 *
 *  Created on: May 7, 2021
 *      Author: phaentu
 */

#ifndef BAM_TYPES_H_
#define BAM_TYPES_H_

#include <cstdint>
#include <ostream>
#include <math.h>
#include <cmath>
#include "strongTypes.h"

namespace BAM{

//------------------------------------------------
// Base
//------------------------------------------------

enum BaseEnum : uint8_t {A=0, C, G, T, N};

class Base : public StrongTypes::StrongType<BaseEnum, Base> {
private:
	static constexpr char _toChar[5] = {'A', 'C', 'G', 'T', 'N'};
	static constexpr BaseEnum _flipBase[5] = {T, G, C, A, N};

	BaseEnum constexpr _toBaseEnum(const char & base){
		if(base == 'A' || base == 'a') return A;
		else if(base == 'C' || base == 'c') return C;
		else if(base == 'G' || base == 'g') return G;
		else if(base == 'T' || base == 't') return T;
		else return N;
	};

public:
	//constructors
	explicit constexpr Base() : StrongType(N) {};
	explicit constexpr Base(const BaseEnum & base) : StrongType(base) {};
	explicit constexpr Base(const char & base){
		_value = _toBaseEnum(base);
	};

	//assignments
	constexpr void operator=(const BaseEnum & base){
		_value = base;
	};

	constexpr void operator=(const char & base){
		_value = _toBaseEnum(base);
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

	constexpr GenotypeEnum _toGenoOrdered(const Base & first, const Base & second){
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

	constexpr GenotypeEnum _toGeno(const Base & first, const Base & second){
		if(static_cast<uint8_t>( (BaseEnum) first) > static_cast<uint8_t>( (BaseEnum) second)){
			return _toGenoOrdered(second, first);
		} else {
			return _toGenoOrdered(first, second);
		}
	};

public:
	//constructors
	explicit constexpr Genotype() : StrongType(NN) {};

	explicit constexpr Genotype(const GenotypeEnum & genotype) : StrongType(genotype) {};

	explicit constexpr Genotype(const Base & first, const Base & second){
		if(static_cast<uint8_t>( (BaseEnum) first) > static_cast<uint8_t>( (BaseEnum) second)){
			_value = _toGenoOrdered(second, first);
		} else {
			_value = _toGenoOrdered(first, second);
		}
	};

	explicit constexpr Genotype(const BaseEnum & first, const BaseEnum & second) : Genotype(Base(first), Base(second)) {};

	//assign
	void set(const BaseEnum & first, const BaseEnum & second){

	};

	//convert
	explicit operator std::string() const {
		return _toString[static_cast<uint8_t>( (GenotypeEnum) _value)];
	};
};

std::ostream& operator<<(std::ostream& os, const Genotype & genotype);

//------------------------------------------------
// BaseContext
//------------------------------------------------

enum BaseContextEnum : uint8_t {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT, cNN}; //N means unknown base or "nothing", i.e. end of read

class BaseContext : public StrongTypes::StrongType<BaseContextEnum, BaseContext> {
private:
	static const std::string _toString[cNN+1];

	constexpr BaseContextEnum _toBaseEnum(const BaseEnum & first, const BaseEnum & second){
		if(second == N){
			return cNN;
		} else {
			return static_cast<BaseContextEnum>( static_cast<uint8_t>(N) * static_cast<uint8_t>( (BaseEnum) first) + static_cast<uint8_t>( (BaseEnum) second));
		}
	};

public:
	//constructors
	explicit constexpr BaseContext() : StrongType(cNN) {};

	explicit constexpr BaseContext(const BaseContextEnum & context) : StrongType(context) {};

	explicit constexpr BaseContext(const Base & first, const Base & second){
		_value = _toBaseEnum((BaseEnum) first, (BaseEnum) second);
	};

	explicit constexpr BaseContext(const BaseEnum & first, const BaseEnum & second){
		_value = _toBaseEnum( first, second);
	};

	//assign
	constexpr void set(const Base & first, const Base & second){
		_value = _toBaseEnum((BaseEnum) first, (BaseEnum) second);
	};

	constexpr void set(const BaseEnum & first, const BaseEnum & second){
		_value = _toBaseEnum(first, second);
	};

	constexpr void set(const Base & first, const BaseEnum & second){
		_value = _toBaseEnum( (BaseEnum) first, second);
	};

	constexpr void set(const BaseEnum & first, const Base & second){
		_value = _toBaseEnum(first, (BaseEnum) second);
	};

	//convert
	explicit operator std::string() const {
		return _toString[static_cast<uint8_t>( (BaseContextEnum) _value)];
	};

	explicit operator uint16_t() const {
		return static_cast<uint8_t>( (BaseContextEnum) _value);
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
class BaseQuality;

class ErrorRate : public StrongTypes::StrongType<double, ErrorRate, StrongTypes::Orderable, StrongTypes::Printable>{
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
	explicit ErrorRate(const BaseQuality & quality);
};

//-------------------------------------
// LogErrorRate
// The log of an error rate
//-------------------------------------
class LogErrorRate : public StrongTypes::StrongType<double, LogErrorRate, StrongTypes::Orderable, StrongTypes::Printable>{
private:
	static const double _phredIntToLogError[256];

	template <typename T>
	constexpr double _phredToLogError(const T & phred){
		return static_cast<double>(phred) * -0.2302585092994046; //log(10) / 10
	};

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
	explicit LogErrorRate(const BaseQuality & quality);
};

//------------------------------------------------
// PhredErrorRate
// phreded error = -10 * log_10(error)
//------------------------------------------------
class PhredErrorRate : public StrongTypes::StrongType<double, PhredErrorRate, StrongTypes::Orderable, StrongTypes::Printable> {
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
	explicit PhredErrorRate(const BaseQuality & quality);
};

//------------------------------------------------
// PhredIntError
// phreded error stored as uint8_t
// only valid within [0,255] and truncated outside
//------------------------------------------------
class PhredIntErrorRate : public StrongTypes::StrongType<uint8_t, PhredIntErrorRate, StrongTypes::Orderable, StrongTypes::Printable> {
private:
	constexpr uint8_t _qualityToPhredInt(const char & quality){
		return quality - 33;
	};

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

	explicit PhredIntErrorRate(const BaseQuality & quality);

	//assignments
	constexpr void operator=(const BaseQuality & quality);
	constexpr void operator=(const uint8_t & phredInt){
		_value = phredInt;
	};
};

//------------------------------------------------
// BaseQuality
// quality = PhredIntError + 33
// only valid within [33,126] (matching chars ! until ~) and truncated outside
//------------------------------------------------
class BaseQuality : public StrongTypes::StrongType<char, BaseQuality, StrongTypes::Orderable, StrongTypes::Printable> {
private:
	static const char _min = 33;
	static const char _max = 126;

	constexpr char _toQualityChar(const char& quality){
		if(quality < _min){
			return _min;
		} else if(quality > _max ){
			return _max;
		} else {
			return quality;
		}
	};

public:
	//constructors
	explicit constexpr BaseQuality() : StrongType() {};
	explicit BaseQuality(const ErrorRate & error);
	explicit BaseQuality(const LogErrorRate & error);
	explicit BaseQuality(const PhredErrorRate & error);

	explicit constexpr BaseQuality(const PhredIntErrorRate & error){
		if( (uint8_t) error > 125){
			_value = 126;
		} else {
			_value = (uint8_t) error + 33;
		}
	};

	explicit constexpr BaseQuality(const char & quality){
		_value = _toQualityChar(quality);
	};

	//assign
	constexpr void operator=(const char & quality){
		_value = _toQualityChar(quality);
	};

	//get range information
	constexpr char min() const { return _min; };
	constexpr char max() const { return _max; };
};

}; //end namespace BAM


#endif /* BAM_TYPES_H_ */
