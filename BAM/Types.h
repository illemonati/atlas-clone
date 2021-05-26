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
#include "stringFunctions.h"
#include "probability.h"
#include "TRandomGenerator.h"

namespace BAM{

//------------------------------------------------
// Base
//------------------------------------------------

enum BaseEnum : uint8_t {A=0, C, G, T, N};

class Base : public StrongTypes::StrongType<BaseEnum, Base> {
private:
	static constexpr char _toChar[5] = {'A', 'C', 'G', 'T', 'N'};
	static constexpr BaseEnum _flipBase[5] = {T, G, C, A, N};

	[[nodiscard]] BaseEnum constexpr _toBaseEnum(const char & base){
		if(base == 'A' || base == 'a') return A;
		else if(base == 'C' || base == 'c') return C;
		else if(base == 'G' || base == 'g') return G;
		else if(base == 'T' || base == 't') return T;
		else return N;
	};

public:
	//constructors
	explicit constexpr Base() : StrongType(N) {};
	constexpr Base(const BaseEnum & base) : StrongType(base) {};
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

	explicit operator std::string() const {
		return std::string( _toChar[ static_cast<uint8_t>(_value) ], 1 );
	};

	explicit constexpr operator uint8_t() const {
		return static_cast<uint8_t>(_value);
	};

	//manipulate
	void flip(){
		_value = _flipBase[ static_cast<uint8_t>(_value) ];
	};
	[[nodiscard]] Base constexpr flipped() const{
		return Base(_flipBase[ static_cast<uint8_t>(_value) ]);
	};

	//range info (to loop)
	[[nodiscard]] static constexpr Base min() { return Base(A); };
	[[nodiscard]] static constexpr Base max() { return Base(N); };
};

std::ostream& operator<<(std::ostream& os, const Base & base);

//------------------------------------------------
// Genotype
//------------------------------------------------

enum GenotypeEnum : uint8_t {AA=0, AC, AG, AT, CC, CG, CT, GG, GT, TT, NN};

class Genotype : public StrongTypes::StrongType<GenotypeEnum, Genotype> {
private:
	static constexpr Base _firstBase[11] = {A, A, A, A, C, C, C, G, G, T, N};
	static constexpr Base _secondBase[11] = {A, C, G, T, C, G, T, G, T, T, N};
	static constexpr bool _isHomozygous[11] = {true, false, false, false, true, false, false, true, false, true, false};
	static constexpr bool _isHeterozygous[11] = {false, true, true, true, false, true, true, false, true, false, false};

	static std::string _toString[NN+1];

	[[nodiscard]] constexpr GenotypeEnum _toGenoOrdered(const Base & first, const Base & second){
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

	[[nodiscard]] constexpr GenotypeEnum _toGeno(const Base & first, const Base & second){
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
		_value = _toGeno(first, second);
	};

	explicit constexpr Genotype(const BaseEnum & first, const BaseEnum & second) : Genotype(Base(first), Base(second)) {};

	//assignments
	constexpr void operator=(const GenotypeEnum & genotype){
		_value = genotype;
	};

	void set(const Base & first, const Base & second){
		_value = _toGeno(first, second);
	};

	void set(const BaseEnum & first, const BaseEnum & second){
		_value = _toGeno(first, second);
	};

	//convert
	[[nodiscard]] explicit operator std::string() const {
		return _toString[static_cast<uint8_t>( (GenotypeEnum) _value)];
	};

	[[nodiscard]] explicit constexpr operator uint8_t() const {
		return static_cast<uint8_t>(_value);
	};

	//get
	[[nodiscard]] constexpr const Base& firstAllele() const{
		return _firstBase[_value];
	};

	[[nodiscard]] constexpr const Base& secondAllele() const{
		return _secondBase[_value];
	};

	[[nodiscard]] constexpr const Base& randomAllele(TRandomGenerator & RandomGenerator) const{
		if(_value == NN){
			return _firstBase[10];
		} else if(_isHomozygous[_value]){
			return firstAllele();
		} else if(RandomGenerator.getRand() < 0.5){
			return firstAllele();
		} else {
			return secondAllele();
		}
	};

	[[nodiscard]] constexpr bool isHomozygous() const{
		return _isHomozygous[_value];
	};

	[[nodiscard]] constexpr bool isHeterozygous() const{
		return _isHeterozygous[_value];
	};

	//++operator and range info (to loop)
	[[nodiscard]] static constexpr Genotype min() { return Genotype(AA); };
	[[nodiscard]] static constexpr Genotype max() { return Genotype(NN); };

	constexpr Genotype& operator++(){
		if(_value == NN){
			throw std::runtime_error("constexpr Genotype& operator++():: overflow!");
		};
		_value = static_cast<GenotypeEnum>( static_cast<uint8_t>(_value) + 1);
		return *this;
	};

	constexpr bool operator<(const Genotype & other){
		return static_cast<uint8_t>(_value) < static_cast<uint8_t>(other.get());
	};
};

std::ostream& operator<<(std::ostream& os, const Genotype & genotype);

//-------------------------------------
// BiallelicGenotypes
//-------------------------------------
class BiallelicGenotypes{
private:
	std::array<Genotype, 3> _genotypes;

public:
	BiallelicGenotypes(const Base & First, const Base & Second){
		_genotypes[0].set(First, First);
		_genotypes[1].set(First, Second);
		_genotypes[2].set(Second, Second);
	};
	~BiallelicGenotypes() = default;

	[[nodiscard]] constexpr const Genotype& homoFirst() const { return _genotypes[0]; };
	[[nodiscard]] constexpr const Genotype& het() const { return _genotypes[1]; };
	[[nodiscard]] constexpr const Genotype& homoSecond() const { return _genotypes[2]; };

	[[nodiscard]] constexpr const Genotype& operator[](const uint8_t& Index) const {
		if(Index > 2){
			throw std::runtime_error("constexpr const Genotype& operator[](const uint8_t& Index) const: Unknown genotype index " + toString(Index) + "!");
		} else {
			return _genotypes[Index];
		}
	};
};

//------------------------------------------------
// BaseContext
//------------------------------------------------

enum BaseContextEnum : uint8_t {cAA=0, cAC, cAG, cAT, cCA, cCC, cCG, cCT, cGA, cGC, cGG, cGT, cTA, cTC, cTG, cTT, cNA, cNC, cNG, cNT, cNN}; //N means unknown base or "nothing", i.e. end of read

class BaseContext : public StrongTypes::StrongType<BaseContextEnum, BaseContext> {
private:
	static const std::string _toString[cNN+1];

	[[nodiscard]] constexpr BaseContextEnum _toBaseContextEnum(const BaseEnum & first, const BaseEnum & second){
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
		_value = _toBaseContextEnum((BaseEnum) first, (BaseEnum) second);
	};

	explicit constexpr BaseContext(const BaseEnum & first, const BaseEnum & second){
		_value = _toBaseContextEnum( first, second);
	};

	//assign
	constexpr void operator=(const BaseContextEnum & context){
		_value = context;
	};

	constexpr void set(const Base & first, const Base & second){
		_value = _toBaseContextEnum((BaseEnum) first, (BaseEnum) second);
	};

	constexpr void set(const BaseEnum & first, const BaseEnum & second){
		_value = _toBaseContextEnum(first, second);
	};

	constexpr void set(const Base & first, const BaseEnum & second){
		_value = _toBaseContextEnum( (BaseEnum) first, second);
	};

	constexpr void set(const BaseEnum & first, const Base & second){
		_value = _toBaseContextEnum(first, (BaseEnum) second);
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

//------------------------------------------------
// AllelicCombination
// stores the combination of alleles considered for genoytpes (e.g. after major / minor)
//------------------------------------------------

enum AllelicCombinationEnum : uint8_t {alleleicCombinationAC=0, alleleicCombinationAG, alleleicCombinationAT, alleleicCombinationCG, alleleicCombinationCT, alleleicCombinationGT, alleleicCombinationNN}; //acNN means unknown

class AllelicCombination : public StrongTypes::StrongType<AllelicCombinationEnum, AllelicCombination> {
private:

	static constexpr Base _firstBase[7] = {A, A, A, C, C, G, N};
	static constexpr Base _secondBase[7] = {C, G, T, G, T, T, N};

	[[nodiscard]] constexpr AllelicCombinationEnum _toAlleleicCombinationEnum(const BaseEnum & first, const BaseEnum & second){
		if(static_cast<uint8_t>( (BaseEnum) first) > static_cast<uint8_t>( (BaseEnum) second)){
			return _toAlleleicCombinationEnum(second, first);
		} else {
			if(first == N || second == N){
				return alleleicCombinationNN;
			} else if(first == second) {
				throw std::runtime_error("constexpr AllelicCombinationEnum _toAlleleicCombinationEnum(const BaseEnum & first, const BaseEnum & second): first == second!");
			} else {
				if(first == A){
					return static_cast<AllelicCombinationEnum>( static_cast<uint8_t>(second) - 1);
				} else if(first ==C){
					return static_cast<AllelicCombinationEnum>( static_cast<uint8_t>(second) + 1);
				} else {
					return alleleicCombinationGT;
				}
			}
		}
	};

public:
	//constructors
	explicit constexpr AllelicCombination() : StrongType(alleleicCombinationNN) {};

	explicit constexpr AllelicCombination(const AllelicCombinationEnum & combination) : StrongType(combination) {};

	explicit constexpr AllelicCombination(const Base & first, const Base & second){
		_value = _toAlleleicCombinationEnum((BaseEnum) first, (BaseEnum) second);
	};

	explicit constexpr AllelicCombination(const BaseEnum & first, const BaseEnum & second){
		_value = _toAlleleicCombinationEnum(first, second);
	};

	//assign
	constexpr void operator=(const AllelicCombinationEnum & combination){
		_value = combination;
	};

	constexpr void set(const Base & first, const Base & second){
		_value = _toAlleleicCombinationEnum((BaseEnum) first, (BaseEnum) second);
	};

	constexpr void set(const BaseEnum & first, const BaseEnum & second){
		_value = _toAlleleicCombinationEnum(first, second);
	};

	constexpr void set(const Base & first, const BaseEnum & second){
		_value = _toAlleleicCombinationEnum( (BaseEnum) first, second);
	};

	constexpr void set(const BaseEnum & first, const Base & second){
		_value = _toAlleleicCombinationEnum(first, (BaseEnum) second);
	};

	//convert / get
	[[nodiscard]] constexpr const Base& firstAllele() const{
		return _firstBase[_value];
	};

	[[nodiscard]] constexpr const Base& secondAllele() const{
		return _secondBase[_value];
	};

	[[nodiscard]] constexpr Genotype homoFirst() const{
		return Genotype(_firstBase[_value], _firstBase[_value]);
	};

	[[nodiscard]] constexpr Genotype het() const{
		return Genotype(_firstBase[_value], _secondBase[_value]);
	};

	[[nodiscard]] constexpr Genotype homoSecond() const{
		return Genotype(_secondBase[_value], _secondBase[_value]);
	};
};

//------------------------------------------------
// PhredErrorRate
// phreded error = -10 * log_10(error)
//------------------------------------------------
class PhredIntErrorRate;
class HighPrecisionPhredIntErrorRate;
class BaseQuality;

class PhredErrorRate : public StrongTypes::StrongType<double, PhredErrorRate, StrongTypes::Orderable, StrongTypes::Printable> {
private:
	[[nodiscard]] constexpr double _errorToPhredError(const double & error) const {
		return -10.0 * log10( (double) error);
	};

	[[nodiscard]] constexpr double _phredErrorToError(const double & phred) const {
		return pow(10.0, phred / -10.0);
	};

	[[nodiscard]] constexpr double _logErrorToPhredError(const double & logError) const {
		return -4.342944819032518 * logError; // -10 * log_10(e) * log(error)
	};

	[[nodiscard]] constexpr double _phredErrorToLogError(const double & logError) const {
		return logError / -4.342944819032518;
	};

	[[nodiscard]] constexpr double _log10ErrorToPhredError(const double & log10Error) const {
		return -10.0 * log10Error;
	};

	[[nodiscard]] constexpr double _phredErrorToLog10Error(const double & log10Error) const {
		return log10Error / -10;
	};

public:
	//constructors
	explicit constexpr PhredErrorRate() : StrongType() {};
	explicit constexpr PhredErrorRate(const double & phredError){
		if(phredError < 0){
			throw std::runtime_error("explicit constexpr PhredErrorRate(const double & phredError): phredError < 0!");
		}
		_value = phredError;
	};

	explicit constexpr PhredErrorRate(const Probability & error){
		_value = _errorToPhredError((double) error);
	};

	explicit constexpr PhredErrorRate(const LogProbability & logError){
		_value = _logErrorToPhredError((double) logError);
	};

	explicit constexpr PhredErrorRate(const Log10Probability & error){
		_value = _log10ErrorToPhredError(error.get());
	};

	explicit PhredErrorRate(const std::string & logError) : StrongType(logError) {};

	explicit PhredErrorRate(const PhredIntErrorRate & error);
	explicit PhredErrorRate(const HighPrecisionPhredIntErrorRate & error);
	explicit PhredErrorRate(const BaseQuality & quality);

	//assign
	constexpr void operator=(const double & phredError){
		_value = phredError;
	};

	constexpr void operator=(const Probability & probability){
		_value = _errorToPhredError( (double) probability);
	};

	constexpr void operator=(const LogProbability & logProbability){
		_value = _logErrorToPhredError( (double) logProbability);
	};

	constexpr void operator=(const Log10Probability & logProbability){
		_value = _log10ErrorToPhredError( (double) logProbability);
	};

	//convert
	explicit operator std::string() const {
		return toString(_value);
	};

	explicit constexpr operator Probability() const {
		return Probability(_phredErrorToError(_value) );
	};

	explicit constexpr operator LogProbability() const {
		return LogProbability( _phredErrorToLogError(_value) );
	};

	explicit constexpr operator Log10Probability() const {
		return Log10Probability( _phredErrorToLog10Error(_value) );
	};

	explicit operator PhredIntErrorRate() const;
	explicit operator HighPrecisionPhredIntErrorRate() const;
	explicit operator BaseQuality() const;
};

//------------------------------------------------
// PhredIntError
// phreded error stored as uint8_t
// only valid within [0,255] and truncated outside
//------------------------------------------------
class PhredIntErrorRate : public StrongTypes::StrongType<uint8_t, PhredIntErrorRate, StrongTypes::Orderable, StrongTypes::Incrementable, StrongTypes::Printable> {
private:
	static const double _phredIntToError[256];

	[[nodiscard]] constexpr uint8_t _qualityToPhredInt(const char & quality) const {
		return quality - 33;
	};

	[[nodiscard]] constexpr uint8_t _phredErrorToPhredInt(const PhredErrorRate & phredError) const {
		if( (double) phredError < 254.5){
			return std::round( (double) phredError);
		} else {
			return 255;
		}
	};

	[[nodiscard]] constexpr double _phredIntToLogError(const double & phred) const {
		return phred / -4.3429448190325181667;
	};

	[[nodiscard]] constexpr double _phredIntToLog10Error(const double & phred) const {
		return phred / -10.0;
	};

public:
	//constructors
	explicit constexpr PhredIntErrorRate() : StrongType() {};
	explicit constexpr PhredIntErrorRate(const uint8_t & phredInt) : StrongType(phredInt) {};

	explicit constexpr PhredIntErrorRate(const Probability & error) : PhredIntErrorRate((PhredErrorRate) error) {};
	explicit constexpr PhredIntErrorRate(const LogProbability & error) : PhredIntErrorRate((PhredErrorRate) error) {};
	explicit constexpr PhredIntErrorRate(const Log10Probability & error) : PhredIntErrorRate((PhredErrorRate) error) {};

	explicit PhredIntErrorRate(const std::string & logError) : StrongType(logError) {};

	explicit constexpr PhredIntErrorRate(const PhredErrorRate & phredError){
		_value = _phredErrorToPhredInt(phredError);
	};

	explicit PhredIntErrorRate(const HighPrecisionPhredIntErrorRate & error);
	explicit PhredIntErrorRate(const BaseQuality & quality);

	//assign
	void operator=(const BaseQuality & quality);

	constexpr void operator=(const uint8_t & phredInt){
		_value = phredInt;
	};

	constexpr void operator=(const Probability & probability){
		_value = _phredErrorToPhredInt(PhredErrorRate(probability));
	};

	constexpr void operator=(const LogProbability & logProbability){
		_value = _phredErrorToPhredInt(PhredErrorRate(logProbability));
	};

	constexpr void operator=(const Log10Probability & logProbability){
		_value = _phredErrorToPhredInt(PhredErrorRate(logProbability));
	};

	constexpr void operator=(const PhredErrorRate & phredError){
		_value = _phredErrorToPhredInt(phredError);
	};

	//convert
	explicit operator std::string() const {
		return toString(_value);
	};

	explicit constexpr operator Probability() const {
		return Probability( _phredIntToError[_value] );
	};

	explicit constexpr operator LogProbability() const {
		return LogProbability( _phredIntToLogError(_value) );
	};

	explicit constexpr operator Log10Probability() const {
		return Log10Probability( _phredIntToLog10Error(_value) );
	};

	explicit operator PhredErrorRate() const;
	explicit operator HighPrecisionPhredIntErrorRate() const;
	explicit operator BaseQuality() const;

	//get range information
	[[nodiscard]] static constexpr PhredIntErrorRate min() { return PhredIntErrorRate(0); };
	[[nodiscard]] static constexpr PhredIntErrorRate max() { return PhredIntErrorRate(255); };

	//manipulate
	void constexpr makeIllumina(){
		if(_value < 35){
			_value = 33;
		} else if(_value < 43){
			_value = 39;
		} else if(_value < 53){
			_value = 48;
		} else if(_value < 58){
			_value = 55;
		} else if(_value < 63){
			_value = 60;
		} else if(_value < 68){
			_value = 66;
		} else if(_value < 72){
			_value = 70;
		} else {
			_value = 73;
		}
	};
};

//------------------------------------------------
// HighPrecisionPhredIntErrorRate
// phred error with more (100 times) precision than PhredIntError
// HighPrecisionPhredIntErrorRate = -1000 * log10(error)
//------------------------------------------------
class HighPrecisionPhredIntErrorRate : public StrongTypes::StrongType<uint16_t, HighPrecisionPhredIntErrorRate, StrongTypes::Orderable, StrongTypes::Printable> {
private:
	static const double _highPrecisionPhredIntToError[65536];

	[[nodiscard]] constexpr uint16_t _phredErrortoHighPrecisionPhredInt(const PhredErrorRate & error) const {
		double tmp = (double) error * 100.0;
		if( tmp < 65534.5){
			return (uint16_t) std::round( tmp );
		} else {
			return 65535;
		}
	};

	[[nodiscard]] constexpr double _highPrecisionPhredIntToLogError(const double & phred) const {
		return phred / -434.29448190325181667;
	};

	[[nodiscard]] constexpr double _highPrecisionPhredIntToLog10Error(const double & phred) const {
		return phred / -1000.0;
	};

public:
	//constructors
	explicit constexpr HighPrecisionPhredIntErrorRate() : StrongType() {};
	explicit constexpr HighPrecisionPhredIntErrorRate(const uint16_t & highPrecisionPhredInt) : StrongType(highPrecisionPhredInt) {};

	explicit constexpr HighPrecisionPhredIntErrorRate(const Probability & error) : HighPrecisionPhredIntErrorRate((PhredErrorRate) error) {};
	explicit constexpr HighPrecisionPhredIntErrorRate(const LogProbability & error) : HighPrecisionPhredIntErrorRate((PhredErrorRate) error) {};
	explicit constexpr HighPrecisionPhredIntErrorRate(const Log10Probability & error) : HighPrecisionPhredIntErrorRate((PhredErrorRate) error) {};

	explicit HighPrecisionPhredIntErrorRate(const std::string & logError) : StrongType(logError) {};

	explicit constexpr HighPrecisionPhredIntErrorRate(const PhredErrorRate & error){
		_value = _phredErrortoHighPrecisionPhredInt(error);
	};

	explicit constexpr HighPrecisionPhredIntErrorRate(const PhredIntErrorRate & error){
		_value = 100 * (uint8_t) error;
	};

	explicit HighPrecisionPhredIntErrorRate(const BaseQuality & quality);

	//assign
	constexpr void operator=(const uint16_t & highPrecisionPhredInt){
		_value = highPrecisionPhredInt;
	};

	constexpr void operator=(const Probability & probability){
		_value = _phredErrortoHighPrecisionPhredInt(PhredErrorRate(probability));
	};

	constexpr void operator=(const LogProbability & logProbability){
		_value = _phredErrortoHighPrecisionPhredInt(PhredErrorRate(logProbability));
	};

	constexpr void operator=(const Log10Probability & log10Probability){
		_value = _phredErrortoHighPrecisionPhredInt(PhredErrorRate(log10Probability));
	};

	constexpr void operator=(const PhredErrorRate & phredError){
		_value = _phredErrortoHighPrecisionPhredInt(phredError);
	};

	//convert
	explicit operator std::string() const {
		return toString(_value);
	};

	explicit constexpr operator Probability() const {
		return Probability( _highPrecisionPhredIntToError[_value] );
	};

	explicit constexpr operator LogProbability() const {
		return LogProbability( _highPrecisionPhredIntToLogError(_value) );
	};

	explicit constexpr operator Log10Probability() const {
		return Log10Probability( _highPrecisionPhredIntToLog10Error(_value) );
	};

	explicit operator PhredErrorRate() const;
	explicit operator PhredIntErrorRate() const;
	explicit operator BaseQuality() const;

	//get range information
	static constexpr HighPrecisionPhredIntErrorRate min() { return HighPrecisionPhredIntErrorRate(0);   };
	static constexpr HighPrecisionPhredIntErrorRate max() { return HighPrecisionPhredIntErrorRate(65535); };

	//manipulate / normalize
	[[nodiscard]] constexpr HighPrecisionPhredIntErrorRate operator-(HighPrecisionPhredIntErrorRate const & other) const {
		if(other > *this){
			throw std::runtime_error("constexpr HighPrecisionPhredIntErrorRate operator-(HighPrecisionPhredIntErrorRate const & other) const: other > *this!");
		}
		return HighPrecisionPhredIntErrorRate(_value - other.get());
	};

	constexpr HighPrecisionPhredIntErrorRate& operator-=(HighPrecisionPhredIntErrorRate const& other){
		if(other > *this){
			throw std::runtime_error("constexpr HighPrecisionPhredIntErrorRate operator-(HighPrecisionPhredIntErrorRate const & other) const: other > *this!");
		}
		_value -= other.get();
		return *this;
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

	[[nodiscard]] constexpr char _toQualityChar(const char& quality) const {
		if(quality < _min){
			return _min;
		} else if(quality > _max ){
			return _max;
		} else {
			return quality;
		}
	};

	[[nodiscard]] constexpr char _toQualityChar(const PhredIntErrorRate& error) const {
		if( (uint8_t) error > 125){
			return 126;
		} else {
			return (uint8_t) error + 33;
		}
	};

public:
	//constructors
	explicit constexpr BaseQuality() : StrongType() {};

	explicit constexpr BaseQuality(const Probability & error) : BaseQuality((PhredIntErrorRate) error) {};
	explicit constexpr BaseQuality(const LogProbability & error) : BaseQuality((PhredIntErrorRate) error) {};
	explicit constexpr BaseQuality(const Log10Probability & error) : BaseQuality((PhredIntErrorRate) error) {};
	explicit constexpr BaseQuality(const PhredErrorRate & error) : BaseQuality( PhredIntErrorRate(error) ) {};

	explicit constexpr BaseQuality(const PhredIntErrorRate & error){
		_value = _toQualityChar(error);
	};

	explicit constexpr BaseQuality(const HighPrecisionPhredIntErrorRate & error){
		_value = _toQualityChar(PhredIntErrorRate(error));
	};

	explicit constexpr BaseQuality(const char & quality){
		_value = _toQualityChar(quality);
	};

	explicit constexpr BaseQuality(const std::string & quality){
		if(quality.size() != 1){
			throw std::runtime_error("explicit constexpr BaseQuality(const std::string & quality): quality is not of size 1!");
		}
		_value = _toQualityChar(quality[0]);
	};

	//assign
	constexpr void operator=(const char & quality){
		_value = _toQualityChar(quality);
	};

	constexpr void operator=(const PhredIntErrorRate & error){
		_value = _toQualityChar(error);
	};

	//convert
	explicit operator std::string() const {
		return std::string(_value, 1);
	};

	explicit constexpr operator Probability() const {
		return (Probability) PhredIntErrorRate(_value);
	};

	explicit constexpr operator LogProbability() const {
		return (LogProbability) PhredIntErrorRate(_value);
	};

	explicit constexpr operator Log10Probability() const {
		return (Log10Probability) PhredIntErrorRate(_value);
	};

	explicit operator PhredErrorRate() const;
	explicit operator PhredIntErrorRate() const;
	explicit operator HighPrecisionPhredIntErrorRate() const;

	//get range information
	static constexpr BaseQuality min() { return BaseQuality(_min); };
	static constexpr BaseQuality max() { return BaseQuality(_max); };

	//manipulate
	void makeIllumina(){
		PhredIntErrorRate tmp(*this);
		tmp.makeIllumina();
		*this = tmp;
	};
};

}; //end namespace BAM


#endif /* BAM_TYPES_H_ */
