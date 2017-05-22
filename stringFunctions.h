/*
 * stringFunctions.h
 *
 *  Created on: May 4, 2012
 *      Author: wegmannd
 */

#ifndef STRINGFUNCTIONS_H_
#define STRINGFUNCTIONS_H_

#include <string>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdio.h>

std::string toString(const int & input);
std::string toString(const long & input);
//std::string toString(const std::vector<int>::size_type & input);
std::string toString(const unsigned int & input);
std::string toString(const unsigned long & input);
std::string toString(const float & input);
std::string toString(const double & input);

int stringToInt(const std::string & s);
long stringToLong(const std::string & s);
double stringToDouble(const std::string & s);
float stringToFloat(const std::string & s);

int stringToIntCheck(const std::string & s);
long stringToLongCheck(const std::string & s);
double stringToDoubleCheck(const std::string & s);
float stringToFloatCheck(const std::string & s);

//check
bool stringContains(const std::string & haystack, std::string needle);
bool stringContainsAny(const std::string & haystack, std::string needle);
bool stringContains(const std::string & haystack, char needle);
bool stringContainsOnly(const std::string & haystack, std::string needle);
bool stringContainsLetters(const std::string & haystack);
bool stringContainsNumbers(const std::string & haystack);
bool allEntriesAreUnique(const std::vector<std::string> vec);
std::string getFirstNonUniqueString(std::vector<std::string> vec);

//modify
void eraseAllOccurences(std::string & s, std::string delim);
void eraseAllOccurencesAny(std::string & s, std::string delim);
void eraseAllWhiteSpaces(std::string & s);

//manipulations
std::string extractBefore(std::string & s, std::string delim);
std::string extractBefore(std::string & s, char delim);
std::string extractBeforeAnyOf(std::string & s, std::string delim);
std::string extractBeforeDoubleSlash(std::string & s);
std::string extractBeforeWhiteSpace(std::string & s);
std::string extractUntil(std::string & s, std::string delim);
std::string extractUntil(std::string & s, char delim);

std::string extractBeforeLast(std::string & s, std::string delim);
std::string extractBeforeLast(std::string & s, char delim);
std::string extractUntilLast(std::string & s, std::string delim);
std::string extractUntilLast(std::string & s, char delim);

std::string extractAfter(std::string & s, std::string delim);
std::string extractAfter(std::string & s, char delim);

std::string extractAfterLast(std::string & s, std::string delim);
std::string extractAfterLast(std::string & s, char delim);

std::string extractPath(std::string & s);

//read
std::string readAfterLast(std::string & s, char delim);

//manipulate
void trimString(std::string & s);
void trimString(std::string & s, std::string what);
void concatenateString(std::vector<std::string> & vec, std::string & s);
void concatenateString(std::vector<std::string> & vec, std::string & s, int from);
void concatenateString(std::vector<std::string> & vec, std::string & s, std::string delim);
void concatenateString(std::vector<int> & vec, std::string & s, std::string delim);
void concatenateString(std::vector<double> & vec, std::string & s, std::string delim);
std::string concatenateString(std::vector<std::string> & vec, std::string delim);
std::string concatenateString(std::vector<int> & vec, std::string delim);
std::string concatenateString(std::vector<double> & vec, std::string delim);

//read
std::string readAfterLast(std::string & s, char delim);

//split into vector
void fillVectorFromString(std::string s, std::vector<std::string> & vec, std::string delim);
void fillVectorFromStringAny(std::string s, std::vector<std::string> & vec, std::string delim);
void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<std::string> & vec, std::string delim);
bool fillVectorFromStringAnySkipEmptyArray(std::string s, double* array, int & size, std::string & delim);
void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<std::string> & vec);
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<std::string> & vec);
bool fillVectorFromStringWhiteSpaceSkipEmptyArray(const std::string & s, double* array, int & size);
void fillVectorFromString(std::string s, std::vector<std::string> & vec, char delim);

void fillVectorFromString(std::string s, std::vector<float> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<double> & vec, char delim);
void fillVectorFromStringAny(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnyCheck(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<double> & vec);

void fillVectorFromString(std::string s, std::vector<int> & vec, char delim);
void fillVectorFromStringAny(std::string s, std::vector<int> & vec, std::string delim);
void fillVectorFromStringAnyCheck(std::string s, std::vector<int> & vec, std::string delim);
void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<int> & vec, std::string delim);
void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<int> & vec, std::string delim);
void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<int> & vec);
void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<int> & vec);
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<int> & vec);
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<int> & vec);

void fillVectorFromString(std::string s, std::vector<long> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<bool> & vec, char delim);
bool fillSequenceFromString(std::string s, std::vector<int> & vec, char delim);

//read from file
void readHeaderAndValues(std::string & filename, std::vector<std::string> & header, std::vector<double> & values);
void readHeaderAndValuesUnique(std::string & filename, std::vector<std::string> & header, std::vector<double> & values);

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, std::string delim);
void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, char delim);
void fillVectorFromLineAny(std::ifstream & is, std::vector<std::string> & vec, std::string delim);
void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<std::string> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<std::string> & vec);

void fillVectorFromLine(std::ifstream & is, std::vector<double> & vec, char delim);
void fillVectorFromLineAny(std::ifstream & is, std::vector<double> & vec, std::string delim);
void fillVectorFromLineAnyCheck(std::ifstream & is, std::vector<double> & vec, std::string delim);
void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<int> & vec);

std::string stringReplace(char needle, std::string replace, std::string & haystack);
std::string stringReplace(std::string needle, std::string replace, std::string haystack);

//expand / repeat index
bool addRepeatedIndexIfRepeated(std::string & orig, std::vector<std::string> & vec);
bool addExpandedIndexIfToExpand(std::string & orig, std::vector<std::string> & vec);
void addRepeatedIndex(std::string & orig, std::vector<std::string> & vec);
void addExpandIndex(std::string & orig, std::vector<std::string> & vec);
void addRepeatedAndExpandIndexes(std::string & orig, std::vector<std::string> & vec);
void repeatAndExpandIndexes(std::vector<std::string> & orig, std::vector<std::string> & vec);
void repeatIndexes(std::vector<std::string> & orig, std::vector<std::string> & vec);
void repeatIndexes(std::vector<std::string> & orig, std::vector<double> & vec);
void repeatIndexes(std::vector<std::string> & orig, std::vector<int> & vec);
void repeatIndexes(std::vector<std::string> & orig, std::vector<long> & vec);
void addRepeatedAndExpandedIndexesOfSub(const std::string & orig, std::vector< std::vector<std::string> > & vec, std::string delim);
void repeatAndExpandIndexesOfSubs(std::vector<std::string> & orig, std::vector< std::vector<std::string> > & vec, std::string delim);

#endif /* STRINGFUNCTIONS_H_ */
