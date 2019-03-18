
#include "stringFunctions.h"

//-----------------------------------------------------------------------
// casting
//-----------------------------------------------------------------------
/*
std::string toString(const int & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const long & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const std::vector<int>::size_type & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const unsigned int & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const unsigned long & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};


std::string toString(const float & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};


std::string toString(const double & input){
	return std::to_string(input);
};
*/

int stringToInt(const std::string & s){
	return atoi(s.c_str());
};

long stringToLong(const std::string & s){
	return atol(s.c_str());
};

double stringToDouble(const std::string & s){
	return atof(s.c_str());
};

float stringToFloat(const std::string & s){
	return atof(s.c_str());
};
int stringToIntCheck(const std::string & s){
	if(!stringContainsOnly(s, "1234567890.Ee-+")) throw "String '" + s + "' is not a number!";
	return (int) stringToLongCheck(s);
};

long stringToLongCheck(const std::string & s){
	if(!stringContainsOnly(s, "1234567890.Ee-+")) throw "String '" + s + "' is not a number!";
	char** ptr=NULL;
	long i=strtol(s.c_str(), ptr, 10);
	if(ptr!=NULL) throw "String '" +s+"' is not a number!";
	return i;
};

double stringToDoubleCheck(const std::string & s){
	if(!stringContainsOnly(s, "1234567890.Ee-+")) throw "String '" + s + "' is not a number!";
	char** ptr=NULL;
	double i=strtod(s.c_str(), ptr);
	if(ptr!=NULL) throw "String '" +s+"' is not a number!";
	return i;
};

float stringToFloatCheck(const std::string & s){
	return (float) stringToDoubleCheck(s);
};


//-----------------------------------------------------------------------
//check
//-----------------------------------------------------------------------
bool stringContains(const std::string & haystack, std::string needle){
	if(haystack.find(needle)==std::string::npos) return false;
	else return true;
};
bool stringContainsAny(const std::string & haystack, std::string needle){
	if(haystack.find_first_of(needle)==std::string::npos) return false;
	else return true;
};
bool stringContains(const std::string & haystack, char needle){
	if(haystack.find_first_of(needle)==std::string::npos) return false;
	else return true;
};
bool stringContainsOnly(const std::string & haystack, std::string needle){
	if(haystack.find_first_not_of(needle)==std::string::npos) return true;
	else return false;
};
bool stringIsProbablyANumber(const std::string & haystack){
	if(stringContainsOnly(haystack, "1234567890.Ee-+")) return true;
	return false;
};
bool stringContainsLetters(const std::string & haystack){
	return stringContainsAny(haystack, "abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVWXYZäöüÄÖÜàéèÀÉÈ");
};
bool stringContainsNumbers(const std::string & haystack){
	return stringContainsAny(haystack, "1234567890");
};
bool allEntriesAreUnique(std::vector<std::string> vec){
	std::vector<std::string>::iterator it_second;
	for(std::vector<std::string>::iterator it=vec.begin(); it!= --vec.end(); ++it){
		it_second=it;
		++it_second;
		for(;it_second!=vec.end(); ++it){
			if((*it).compare(*it_second)==0) return false;
		}
	}
	return true;
}
std::string getFirstNonUniqueString(std::vector<std::string> vec){
	std::vector<std::string>::iterator it_second;
	for(std::vector<std::string>::iterator it=vec.begin(); it!= --vec.end(); ++it){
		it_second=it;
		++it_second;
		for(;it_second!=vec.end(); ++it_second){
			if((*it).compare(*it_second)==0) return *it;
		}
	}
	return "";
}

//-----------------------------------------------------------------------
//modify
//-----------------------------------------------------------------------
void eraseAllOccurences(std::string & s, std::string delim){
	std::string::size_type l=s.find(delim);
	while(l!=std::string::npos){
		s.erase(l, 1);
		l=s.find_first_of(delim);
	}
};

void eraseAllOccurencesAny(std::string & s, std::string delim){
	std::string::size_type l=s.find_first_of(delim);
	while(l!=std::string::npos){
		s.erase(l, 1);
		l=s.find_first_of(delim);
	}
};

void eraseAllWhiteSpaces(std::string & s){
	eraseAllOccurencesAny(s, " \t\f\v\n\r");
};

//-----------------------------------------------------------------------
//extract before
//-----------------------------------------------------------------------
std::string extractBefore(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBefore(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find_first_of(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeAnyOf(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find_first_of(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeDoubleSlash(std::string & s){
	return extractBefore(s, "//");
};

std::string extractBeforeWhiteSpace(std::string & s){
	return extractBeforeAnyOf(s, " \t\f\v\n\r");
};

std::string extractUntil(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntil(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeLast(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntilLast(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntilLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

//-----------------------------------------------------------------------
//extract after
//-----------------------------------------------------------------------
std::string extractAfter(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfter(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfterLast(std::string & s, std::string delim) {
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfterLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractPath(std::string & s){
	return extractUntilLast(s, "/");
};

//-----------------------------------------------------------------------
//read
//-----------------------------------------------------------------------
std::string readAfterLast(const std::string & s, char delim){
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos) return s.substr(l+1);
	else return s;
};

std::string readBeforeLast(const std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
	} else {
		ret=s;
	}
	return ret;
};

std::string readBeforeLast(const std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
	} else {
		ret=s;
	}
	return ret;
};

//-----------------------------------------------------------------------
//manipulate
//-----------------------------------------------------------------------
void trimString(std::string & s){
	trimString(s, " \t\f\v\n\r");
};

void trimString(std::string & s, std::string what){
	//from beginning
	std::string::size_type l=s.find_first_not_of(what);
	if(l==std::string::npos){
		s.clear();
	} else {
		s.erase(0, l);
		//from end
		l=s.find_last_not_of(what);
		if(l!=std::string::npos)
			s.erase(l+1);
	}
};

void concatenateString(std::vector<std::string> & vec, std::string & s){
	s.clear();
	if(vec.size()>0){
		for(std::vector<std::string>::iterator it=vec.begin(); it!=vec.end(); ++it){
			s+=(*it);
		}
	}
};
void concatenateString(std::vector<std::string> & vec, std::string & s, int from){
	s.clear();
	if(vec.size()>0){
		std::vector<std::string>::iterator it=vec.begin();
		it+=from;
		for(; it!=vec.end(); ++it){
			s+=(*it);
		}
	}
}
void concatenateString(std::vector<std::string> & vec, std::string & s, std::string delim){
	s.clear();
	if(vec.size()>0){
		std::vector<std::string>::iterator it=vec.begin();
		s=*it;
		++it;
		for(; it!=vec.end(); ++it){
			s+=delim+(*it);
		}
	}
};

void concatenateString(std::vector<int> & vec, std::string & s, std::string delim){
	s.clear();
	if(vec.size()>0){
		std::vector<int>::iterator it=vec.begin();
		s=toString(*it);
		++it;
		for(; it!=vec.end(); ++it){
			s+=delim+toString(*it);
		}
	}
}

void concatenateString(std::vector<double> & vec, std::string & s, std::string delim){
	s.clear();
	if(vec.size()>0){
		std::vector<double>::iterator it=vec.begin();
		s=toString(*it);
		++it;
		for(; it!=vec.end(); ++it){
			s+=delim+toString(*it);
		}
	}
}
std::string concatenateString(std::vector<std::string> & vec, std::string delim){
	std::string s;
	concatenateString(vec, s, delim);
	return s;
}
std::string concatenateString(std::vector<int> & vec, std::string delim){
	std::string s;
	if(vec.size()>0){
		std::vector<int>::iterator it=vec.begin();
		s=toString(*it);
		++it;
		for(; it!=vec.end(); ++it){
			s+=delim+toString(*it);
		}
	}
	return s;
}
std::string concatenateString(std::vector<double> & vec, std::string delim){
	std::string s;
	if(vec.size()>0){
		std::vector<double>::iterator it=vec.begin();
		s=toString(*it);
		++it;
		for(; it!=vec.end(); ++it){
			s+=delim+toString(*it);
		}
	}
	return s;
}

std::string concatenateString(double* array, int length, std::string delim){
	std::string s;
	if(length>0){
		s = toString(array[0]);
		for(int i=1; i<length; ++i){
			s += delim + toString(array[i]);
		}
	}
	return s;
}

//-----------------------------------------------------------------------
//split into vector
//-----------------------------------------------------------------------
void fillVectorFromString(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromStringAny(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(tmp);
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty()) vec.push_back(tmp);
	}
};

void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<std::string> & vec){
	fillVectorFromStringAny(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<std::string> & vec){
	fillVectorFromStringAnySkipEmpty(s, vec, " \t\f\v\n\r");
}

//--------------------------------------------------------------------------------------------------------------------
//fillVectorFromString
//--------------------------------------------------------------------------------------------------------------------
void fillVectorFromString(std::string s, std::vector<std::string> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromString(std::string s, std::vector<int> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atoi(s.substr(0,l).c_str()));
	}
};

void fillVectorFromString(std::string s, std::vector<long> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToLong(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atol(s.substr(0,l).c_str()));
	}
};

void fillVectorFromString(std::string s, std::vector<bool> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atoi(s.substr(0,l).c_str()));
	}
};

void fillVectorFromString(std::string s, std::vector<float> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromString(std::string s, std::vector<double> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

//--------------------------------------------------------------------------------------------------------------------
//fillVectorFromStringSkipEmpty
//--------------------------------------------------------------------------------------------------------------------
void fillVectorFromStringSkipEmpty(std::string s, std::vector<std::string> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(tmp);
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty()) vec.push_back(tmp);
	}
};

void fillVectorFromStringSkipEmpty(std::string s, std::vector<int> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty())
				vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty())
			vec.push_back(stringToInt(s.substr(0,l)));
	}
};

void fillVectorFromStringSkipEmpty(std::string s, std::vector<long> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty())
				vec.push_back(stringToLong(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty())
			vec.push_back(stringToLong(s.substr(0,l)));
	}
};

void fillVectorFromStringSkipEmpty(std::string s, std::vector<bool> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty())
				vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty())
			vec.push_back(stringToInt(s.substr(0,l)));
	}
};

void fillVectorFromStringSkipEmpty(std::string s, std::vector<float> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty())
				vec.push_back(stringToFloat(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty())
			vec.push_back(stringToFloat(s.substr(0,l)));
	}
};

void fillVectorFromStringSkipEmpty(std::string s, std::vector<double> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp = s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty())
				vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp = s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty())
			vec.push_back(stringToDouble(s.substr(0,l)));
	}
};

//--------------------------------------------------------------------------------------------------------------------
//fillVectorFromStringAny
//--------------------------------------------------------------------------------------------------------------------
void fillVectorFromStringAny(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromStringAnyCheck(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDoubleCheck(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDoubleCheck(s.substr(0)));
	}
};

void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToDouble(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToDouble(s.substr(0)));
	}
};
bool fillVectorFromStringAnySkipEmptyArray(std::string s, double* array, int & size, std::string delim){
	if(s.empty()) return false;
	std::string::size_type l=s.find_first_of(delim);
	std::string tmp;
	int i = 0;
	while(l!=std::string::npos){
		tmp=s.substr(0,l);
		trimString(tmp);
		if(!tmp.empty()){
			if(i == size) return false;
			array[i] = stringToDouble(tmp);
			++i;
		}
		s.erase(0, l+1);
		l=s.find_first_of(delim);
	}
	tmp=s.substr(0);
	if(!tmp.empty()){
		if(i == size) return false;
		array[i] = stringToDouble(tmp);
		++i;
	}
	if(i != size) return false;
	return true;
};
void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToDouble(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToDoubleCheck(s.substr(0)));
	}
};

void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAny(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnyCheck(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnySkipEmpty(s, vec, " \t\f\v\n\r");
};
bool fillVectorFromStringWhiteSpaceSkipEmptyArray(const std::string & s, double* array, int & size){
	return fillVectorFromStringAnySkipEmptyArray(s, array, size, " \t\f\v\n\r");
}
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnySkipEmptyCheck(s, vec, " \t\f\v\n\r");
};

void fillVectorFromStringAny(std::string s, std::vector<int> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);		}
		vec.push_back(stringToInt(s.substr(0)));
	}
}

void fillVectorFromStringAnyCheck(std::string s, std::vector<int> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToIntCheck(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);		}
		vec.push_back(stringToIntCheck(s.substr(0)));
	}
}
void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<int> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToInt(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToInt(s.substr(0)));
	}
}
void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<int> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToIntCheck(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToIntCheck(s.substr(0)));
	}
}

void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<int> & vec){
	fillVectorFromStringAny(s, vec, " \t\f\v\n\r");
}

void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<int> & vec){
	fillVectorFromStringAnyCheck(s, vec, " \t\f\v\n\r");
}

void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<int> & vec){
	fillVectorFromStringAnySkipEmpty(s, vec, " \t\f\v\n\r");
}
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<int> & vec){
	fillVectorFromStringAnySkipEmptyCheck(s, vec, " \t\f\v\n\r");
}

bool fillSequenceFromString(std::string s, std::vector<int> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::vector<std::string> temp;
		fillVectorFromString(s, temp, delim);
		for(std::vector<std::string>::iterator it=temp.begin(); it!=temp.end(); ++it){
		   //If sequence fill sequence...
			std::string::size_type pos=it->find_first_of('-');
		   if(pos != std::string::npos){
			   int first=atoi((it->substr(0,pos).c_str()));
			   int second=atoi((it->substr(pos+1).c_str()));
			   if(second>first){
				   for(int j=first; j<=second; ++j) vec.push_back(j);
			   } else return false;
		   }
		   //if number, put back.
		   else {
			   int num=atoi(it->c_str());
			   vec.push_back(num);
		   }
		}
		return true;
	} else return false;
};

//-----------------------------------------------------------------------
//read from file
//-----------------------------------------------------------------------
void readHeaderAndValues(std::string & filename, std::vector<std::string> & header, std::vector<double> & values){
	std::string line;
	header.clear();
	values.clear();
	//open file stream
	std::ifstream is (filename.c_str()); // opening the file for reading
	if(!is) throw "The file '" + filename + "' could not be opened!";

	//read header
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, header);

	//read observed Data
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, values);

	//close stream
	is.close();

	if(header.size() != values.size())
		throw "Number of values does not match number of header names in file '" + filename + "'!";
}
void readHeaderAndValuesUnique(std::string & filename, std::vector<std::string> & header, std::vector<double> & values){
	std::string line;
	header.clear();
	values.clear();
	//open file stream
	std::ifstream is (filename.c_str()); // opening the file for reading
	if(!is) throw "The file '" + filename + "' could not be opened!";

	//read header
	getline(is, line);
	trimString(line);
	if(line.empty()) throw "The file '"+filename+"' appears to be empty!";
	fillVectorFromStringWhiteSpaceSkipEmpty(line, header);
	std::string s=getFirstNonUniqueString(header);
	if(!s.empty()) throw "Entry '"+s+"' is listed multiple times in header of file '"+ filename +"'!";

	//read observed Data
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, values);

	//close stream
	is.close();

	if(header.size() != values.size())
		throw "Number of values does not match number of header names in file '" + filename + "'!";
}

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, char delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLineAny(std::ifstream & is, std::vector<std::string> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAny(line, vec, delim);
};

void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<std::string> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpace(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<std::string> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
};

void fillVectorFromLine(std::ifstream & is, std::vector<double> & vec, char delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLineAny(std::ifstream & is, std::vector<double> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAny(line, vec, delim);
};
void fillVectorFromLineAnyCheck(std::ifstream & is, std::vector<double> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAnyCheck(line, vec, delim);
}

void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpace(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<int> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, vec);
}
//-----------------------------------------------------------------------
//search and replace
//-----------------------------------------------------------------------
std::string stringReplace(char needle, std::string replace, std::string & haystack){
	std::string s="";
	std::string::size_type l=haystack.find_first_of(needle);
	while(l!=std::string::npos){
		s=s+haystack.substr(0,l)+replace;
		haystack.erase(0, l+1);
		l=haystack.find_first_of(needle);
	}
	s=s+haystack;
	return s;
};

std::string stringReplace(std::string needle, std::string replace, std::string haystack){
	std::string s="";
	std::string::size_type l=haystack.find(needle);
	while(l!=std::string::npos){
		s=s+haystack.substr(0,l)+replace;
		haystack.erase(0, l+needle.size());
		l=haystack.find(needle);
	}
	s=s+haystack;
	return s;
};

//-----------------------------------------------------------------------
//Repeat / expand indexes in vector notation
//-----------------------------------------------------------------------
//replace entry blah{x} with x times a blah entry
//replace entry blah[x] with blah_1, ..., blah_x

bool addRepeatedIndexIfRepeated(std::string & orig, std::vector<std::string> & vec){
	std::string::size_type pos = orig.find_last_of('{');
	if(pos != std::string::npos){
		if(orig.find_first_of('{') != pos) throw "Multiple '{' characters in string to repeat '" + orig + "'!";
		if(orig.find_last_of('}') != orig.size()-1)  throw "String to repeat '" + orig + "' does not end with '}'!";
		if(orig.find_first_of('[') != std::string::npos)  throw "String to repeat '" + orig + "' contains a conflicting '[' character!";
		if(orig.find_first_of(']') != std::string::npos)  throw "String to repeat '" + orig + "' contains a conflicting ']' character!";
		std::string tmp = orig.substr(0, pos);
		int len = stringToInt(orig.substr(pos + 1, orig.size() - pos - 2));
		if(len <= 0) throw "Request to repeat string '" + orig + "' zero times!";
		for(int i=1; i<=len; ++i) vec.push_back(tmp);
		return true;
	} else return false;
}

bool addExpandedIndexIfToExpand(std::string & orig, std::vector<std::string> & vec){
	std::string::size_type pos = orig.find_last_of('[');
	if(pos != std::string::npos){
		if(orig.find_first_of('[') != pos) throw "Multiple '[' characters in string to expand '" + orig + "'!";
		std::string::size_type pos2 = orig.find_last_of(']');
		if(pos2 == std::string::npos) throw "Missing closing ']' in string to expand '" + orig + "'!";
		if(orig.find_first_of(']') != pos2) throw "Multiple ']' characters in string to expand '" + orig + "'!";
		if(pos2 < pos) throw "Unable to understand string to expand '" + orig + "': wrong order of '[' and ']'!";
		int len = stringToInt(orig.substr(pos + 1, pos2 - pos));
		if(len <= 0) throw "Request to expand string '" + orig + "' zero times!";
		std::string tmp = orig.substr(0, pos);
		std::string tmp2 = orig.substr(pos2+1);
		for(int i=1; i<=len; ++i) vec.push_back(tmp + toString(i) + tmp2);
		return true;
	} else return false;
}

void addRepeatedIndex(std::string & orig, std::vector<std::string> & vec){
	if(!addRepeatedIndexIfRepeated(orig, vec)) vec.push_back(orig);
}

void addExpandIndex(std::string & orig, std::vector<std::string> & vec){
	if(!addExpandedIndexIfToExpand(orig, vec)) vec.push_back(orig);
}

void addRepeatedAndExpandIndexes(std::string & orig, std::vector<std::string> & vec){
	if(!addRepeatedIndexIfRepeated(orig, vec)){
		if(!addExpandedIndexIfToExpand(orig, vec)) vec.push_back(orig);
	}
}

void repeatAndExpandIndexes(std::vector<std::string> & orig, std::vector<std::string> & vec){
	vec.clear();
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedAndExpandIndexes(*it, vec);
	}
}

void repeatIndexes(std::vector<std::string> & orig, std::vector<std::string> & vec){
	vec.clear();
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedIndex(*it, vec);
	}
}
void repeatIndexes(std::vector<std::string> & orig, std::vector<double> & vec){
	std::vector<std::string> tmp;
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedIndex(*it, tmp);
	}
	for(std::vector<std::string>::iterator it=tmp.begin(); it!=tmp.end(); ++it){
		vec.push_back(stringToDoubleCheck(*it));
	}
}
void repeatIndexes(std::vector<std::string> & orig, std::vector<int> & vec){
	std::vector<std::string> tmp;
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedIndex(*it, tmp);
	}
	for(std::vector<std::string>::iterator it=tmp.begin(); it!=tmp.end(); ++it){
		vec.push_back(stringToIntCheck(*it));
	}
}
void repeatIndexes(std::vector<std::string> & orig, std::vector<long> & vec){
	std::vector<std::string> tmp;
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedIndex(*it, tmp);
	}
	for(std::vector<std::string>::iterator it=tmp.begin(); it!=tmp.end(); ++it){
		vec.push_back(stringToLongCheck(*it));
	}
}

void addRepeatedAndExpandedIndexesOfSub(const std::string & orig, std::vector< std::vector<std::string> > & vec, std::string delim){
	std::vector<std::string> origVec;
	fillVectorFromStringAnySkipEmpty(orig, origVec, delim);
	std::vector<std::string>* tmpVec = new std::vector<std::string>[orig.size()];
	unsigned int times = 1;

	//expand individually
	unsigned int i=0;
	for(std::vector<std::string>::iterator it=origVec.begin(); it!=origVec.end(); ++it, ++i){
		addRepeatedAndExpandIndexes(*it, tmpVec[i]);
		if(tmpVec[i].size() > 1){
			if(times > 1){
				if(tmpVec[i].size()!=times) throw "Unequal number of expansions / repeats in '" + orig + "'!";
			} else times =  tmpVec[i].size();
		}
	}

	//construct new vectors
	for(i=0; i<times; ++i) vec.push_back(std::vector<std::string>());
	for(i=0; i<origVec.size(); ++i){
		if(tmpVec[i].size()==1){
			std::vector< std::vector<std::string> >::reverse_iterator it=vec.rbegin();
			for(unsigned int j=0; j<times; ++j, ++it)
				it->push_back(*(tmpVec[i].begin()));
		} else {
			std::vector< std::vector<std::string> >::reverse_iterator it=vec.rbegin();
			for(std::vector<std::string>::reverse_iterator sIt=tmpVec[i].rbegin(); sIt!=tmpVec[i].rend(); ++sIt, ++it){
				it->push_back(*sIt);
			}
		}
	}
	delete[] tmpVec;
}

void repeatAndExpandIndexesOfSubs(std::vector<std::string> & orig, std::vector< std::vector<std::string> > & vec, std::string delim){
	vec.clear();
	for(std::vector<std::string>::iterator it=orig.begin(); it!=orig.end(); ++it){
		addRepeatedAndExpandedIndexesOfSub(*it, vec, delim);
	}
}

