/*
 * Types.cpp
 *
 *  Created on: May 8, 2021
 *      Author: phaentu
 */

#include "Types.h"


namespace BAM{

//------------------------------------------------
// Base
//------------------------------------------------

std::ostream& operator<<(std::ostream& os, const Base & base){
	if(base == A){
		os << 'A';
	} else if(base == C){
		os << 'C';
	} else if(base == G){
		os << 'G';
	} else if(base == T){
		os << 'T';
	} else {
		os << 'N';
	}
	return os;
};

//------------------------------------------------
// Genotype
//------------------------------------------------
std::string Genotype::_toString[NN+1] = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "NN"};

std::ostream& operator<<(std::ostream& os, const Genotype & genotype){
	os << (std::string) genotype;
	return os;
};

//------------------------------------------------
// Context
//------------------------------------------------
const std::string BaseContext::_toString[cNN+1] = {"AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "NA", "NC", "NG", "NT", "NN"};

std::ostream& operator<<(std::ostream& os, const BaseContext & context){
	os << (std::string) context;
	return os;
};

//------------------------------------------------
// AllelicCombination
//------------------------------------------------
/*
static const Base AllelicCombination::_firstBase[7] = {A, A, A, C, C, G, N};
static const Base AllelicCombination::_secondBase[7] = {C, G, T, G, T, T, N};
*/

//-------------------------------------
// ErrorRate
// An error rate within [0,1]
//-------------------------------------
ErrorRate::ErrorRate(const PhredErrorRate & error){
	_value = std::pow(10.0, - (double) error / 10.0);
};

ErrorRate::ErrorRate(const LogErrorRate & error){
	_value = exp( (double) error);
};

ErrorRate::ErrorRate(const PhredIntErrorRate & error){
	_value = _phredIntToError[ (uint8_t) error];
};

ErrorRate::ErrorRate(const HighPrecisionPhredIntErrorRate & error){
	//exception: use casting operator as translation table is owned by HighPrecisionPhredIntErrorRate
	_value = _highPrecisionPhredIntToError[ (uint16_t) error];
};

ErrorRate::ErrorRate(const BaseQuality & quality){
	_value = _phredIntToError[ (char) quality - 33];
};

ErrorRate::operator LogErrorRate() const {
	return LogErrorRate(*this);
};

ErrorRate::operator PhredErrorRate() const {
	return PhredErrorRate(*this);
};

ErrorRate::operator PhredIntErrorRate() const {
	return PhredIntErrorRate(*this);
};

ErrorRate::operator HighPrecisionPhredIntErrorRate() const{
	return HighPrecisionPhredIntErrorRate(*this);
};

ErrorRate::operator BaseQuality() const{
	return BaseQuality(*this);
};

//-------------------------------------
// LogErrorRate
// The log of an error rate
//-------------------------------------
LogErrorRate::LogErrorRate(const PhredErrorRate & error){
	_value = _phredToLogError((double) error);
};

LogErrorRate::LogErrorRate(const PhredIntErrorRate & error){
	_value = _phredToLogError((uint8_t) error);
};

LogErrorRate::LogErrorRate(const HighPrecisionPhredIntErrorRate & error){
	_value = _phredToLogError((uint16_t) error / 100.0);
};

LogErrorRate::LogErrorRate(const BaseQuality & quality){
	_value = _phredToLogError((char) quality - 33);
};

LogErrorRate::operator ErrorRate() const {
	return ErrorRate(*this);
};

LogErrorRate::operator PhredErrorRate() const {
	return PhredErrorRate(*this);
};

LogErrorRate::operator PhredIntErrorRate() const {
	return PhredIntErrorRate(*this);
};

LogErrorRate::operator HighPrecisionPhredIntErrorRate() const{
	return HighPrecisionPhredIntErrorRate(*this);
};

LogErrorRate::operator BaseQuality() const{
	return BaseQuality(*this);
};

//------------------------------------------------
// PhredErrorRate
// phreded error = -10 * log_10(error)
//------------------------------------------------
PhredErrorRate::PhredErrorRate(const BaseQuality & quality) : PhredErrorRate( PhredIntErrorRate(quality)) {};

PhredErrorRate::PhredErrorRate(const PhredIntErrorRate & error){
	_value = (uint8_t) error;
};

PhredErrorRate::PhredErrorRate(const HighPrecisionPhredIntErrorRate & error){
	_value = (uint16_t) error / 100.0;
};

PhredErrorRate::operator ErrorRate() const {
	return ErrorRate(*this);
};

PhredErrorRate::operator LogErrorRate() const {
	return LogErrorRate(*this);
};

PhredErrorRate::operator PhredIntErrorRate() const {
	return PhredIntErrorRate(*this);
};

PhredErrorRate::operator HighPrecisionPhredIntErrorRate() const{
	return HighPrecisionPhredIntErrorRate(*this);
};

PhredErrorRate::operator BaseQuality() const{
	return BaseQuality(*this);
};

//------------------------------------------------
// PhredIntError
// phreded error stored as uint8_t
// only valid within [0,255] and truncated outside
//------------------------------------------------
PhredIntErrorRate::PhredIntErrorRate(const HighPrecisionPhredIntErrorRate & error){
	if((uint16_t) error > 25500){ //255*100
		_value = 255;
	} else {
		_value = round((uint16_t) error / 100.0);
	}
};

PhredIntErrorRate::PhredIntErrorRate(const BaseQuality & quality){
	_value = _qualityToPhredInt((char) quality);
};

void PhredIntErrorRate::operator=(const BaseQuality & quality){
	_value = _qualityToPhredInt((char) quality);
};

PhredIntErrorRate::operator ErrorRate() const {
	return ErrorRate(*this);
};

PhredIntErrorRate::operator LogErrorRate() const {
	return LogErrorRate(*this);
};

PhredIntErrorRate::operator PhredErrorRate() const {
	return PhredErrorRate(*this);
};

PhredIntErrorRate::operator HighPrecisionPhredIntErrorRate() const{
	return HighPrecisionPhredIntErrorRate(*this);
};

PhredIntErrorRate::operator BaseQuality() const{
	return BaseQuality(*this);
};

//------------------------------------------------
// HighPrecisionPhredIntErrorRate
// phred error with more (100 times) precision than PhredIntError
// HighPrecisionPhredIntErrorRate = -1000 * log10(error)
//------------------------------------------------

HighPrecisionPhredIntErrorRate::HighPrecisionPhredIntErrorRate(const BaseQuality & quality){
	_value = 100 * ((char) quality - 33);
};

HighPrecisionPhredIntErrorRate::operator ErrorRate() const {
	return ErrorRate(*this);
};

HighPrecisionPhredIntErrorRate::operator LogErrorRate() const {
	return LogErrorRate(*this);
};

HighPrecisionPhredIntErrorRate::operator PhredErrorRate() const {
	return PhredErrorRate(*this);
};

HighPrecisionPhredIntErrorRate::operator PhredIntErrorRate() const {
	return PhredIntErrorRate(*this);
};

HighPrecisionPhredIntErrorRate::operator BaseQuality() const {
	return BaseQuality(*this);
};

//------------------------------------------------
// Quality
// Sequencing or mapping quality = PhredIntError + 33
// onyl valid within [33,255] and truncated outside
//------------------------------------------------
BaseQuality::operator ErrorRate() const {
	return ErrorRate(*this);
};

BaseQuality::operator LogErrorRate() const {
	return LogErrorRate(*this);
};

BaseQuality::operator PhredErrorRate() const {
	return PhredErrorRate(*this);
};

BaseQuality::operator PhredIntErrorRate() const {
	return PhredIntErrorRate(*this);
};

BaseQuality::operator HighPrecisionPhredIntErrorRate() const {
	return HighPrecisionPhredIntErrorRate(*this);
};

}; //end namespace BAM
