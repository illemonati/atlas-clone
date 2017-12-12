//---------------------------------------------------------------------------
#include "TParameters.h"
#include <sstream>
//---------------------------------------------------------------------------
TParameters::TParameters(){
	inputFileRead = false;
}

TParameters::TParameters(std::vector<std::string> & commandLineParams, TLog* logfile){
	inputFileRead = false;
	initialize(commandLineParams, logfile);
}
TParameters::TParameters(int & argc, char** argv, TLog* logfile){
	inputFileRead = false;
	std::vector<std::string> commandLineParams;

	//skip first: it is name of executable
	for(int i=1;i<argc;++i)  commandLineParams.push_back(argv[i]);
	initialize(commandLineParams, logfile);
}

void TParameters::addParameter(std::string name, std::string value){
	//check if parameter already exists
	mapParameter[name]  = value;
	parameterUsed[name] = false;
}

void TParameters::initialize(std::vector<std::string> & commandLineParams, TLog* logfile){
	std::string my_name;
	//check if first is name of an input file which means no '='!
	if(!commandLineParams.empty() && !stringContains(commandLineParams[0], '=')){
		readInputfile(commandLineParams[0], logfile);
	}

	//parse command line params and overwrite input file
	std::vector<std::string>::iterator it=commandLineParams.begin();
	if(inputFileRead) ++it;
	for(; it!=commandLineParams.end(); ++it){
		my_name=extractBefore(*it, '=');
		if(stringContains(*it, '=')) mapParameter[my_name]=extractAfter(*it, '=');
		else mapParameter[my_name]="";
	}

	//prepare map to store if a parameter was used
	curParameter=mapParameter.begin();
	endParameter=mapParameter.end();
	for(;curParameter!=endParameter;++curParameter){
		parameterUsed[curParameter->first]=false;
	}
	endParameter=mapParameter.end();
}


//---------------------------------------------------------------------------
void TParameters::readInputfile(std::string fileName, TLog* logfile){
	inputFileName=fileName;
	logfile->listFlush("Reading inputfile '" + (std::string) inputFileName + "' ...");
	std::ifstream is (fileName.c_str());
	if(!is){
		logfile->write(" file could not be opened!");
		logfile->conclude("interpreting '" + (std::string) inputFileName + "' as argument");
	} else {
		std::string buf, my_name, my_value;
		std::string line;

		while(is.good() && !is.eof()){
			getline(is, line);
			line=extractBeforeDoubleSlash(line);
			trimString(line);
			if(!line.empty()){
				my_name=extractBeforeWhiteSpace(line);
				trimString(line);
				my_value=extractBeforeDoubleSlash(line);
				if(!my_name.empty()){
					mapParameter[my_name] = my_value;
				}
			}
		}
		inputFileRead=true;
		logfile->done();
	}
}
//---------------------------------------------------------------------------
bool TParameters::parameterExists(std::string my_name){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return true;
		}
	}
    return false;
}
//---------------------------------------------------------------------------
std::string TParameters::getParameter(std::string & my_name, bool mandatory){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return curParameter->second.c_str();
		}
	}
	if(mandatory)
		throw "The parameter '" + my_name + "' is not defined in the inputfile!";
   else return "";
}
//---------------------------------------------------------------------------
std::string TParameters::getParameterString(std::string my_name, bool mandatory){
	curParameter=mapParameter.begin();
	for(;curParameter!=endParameter;++curParameter){
		if(curParameter->first==my_name){
			parameterUsed[curParameter->first]=true;
			return curParameter->second.c_str();
		}
	}
	if(mandatory){
		if(inputFileRead) throw "The parameter '" + my_name + "' is not defined on the command line nor in the input file '" + inputFileName + "'! ";
		else throw "The parameter '" + my_name + "' is not defined! ";
	}
   else return "";
}
std::string TParameters::getParameterStringWithDefault(std::string my_name, std::string def){
	std::string str=getParameter(my_name, false);
	if(str.empty()) return def;
	else return str;
}
//---------------------------------------------------------------------------
double TParameters::getParameterDouble(std::string my_name, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	if(str.empty()) return 0.0;
	else return stringToDouble(str);
}
double TParameters::getParameterDoubleWithDefault(std::string my_name, double def){
	std::string str=getParameter(my_name, false);
	if(str.empty()) return def;
	else return stringToDouble(str);
}
//---------------------------------------------------------------------------
int TParameters::getParameterInt(std::string my_name, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	if(str.empty()) return 0;
	else return stringToInt(str);
}
int TParameters::getParameterIntWithDefault(std::string my_name, int def){
	std::string str = getParameter(my_name, false);
	if(str.empty()) return def;
	else return stringToInt(str);
}
//---------------------------------------------------------------------------
long TParameters::getParameterLong(std::string my_name, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	if(str.empty()) return 0L;
	else return stringToLong(str);
}
long TParameters::getParameterLongWithDefault(std::string my_name, long def){
	std::string str = getParameter(my_name, false);
	if(str.empty()) return def;
	else return stringToLong(str);
}
//---------------------------------------------------------------------------
//WHY DOES TEMPLATE NOT WORK???
//---------------------------------------------------------------------------
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<int> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<long> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<double> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<float> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<bool> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVector(std::string my_name, std::vector<std::string> & vec, char delim, bool mandatory){
	std::string str=getParameter(my_name, mandatory);
	fillVectorFromString(str, vec, delim);
}
//---------------------------------------------------------------------------
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<int> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<long> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<double> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<float> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<bool> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}
void TParameters::fillParameterIntoVectorWithDefault(std::string my_name, std::vector<std::string> & vec, char delim, std::string def){
	std::string str = getParameterStringWithDefault(my_name, def);
	fillVectorFromString(str, vec, delim);
}

//---------------------------------------------------------------------------
std::string TParameters::getListOfUnusedParameters(){
	std::string parameterList="";
	std::map<std::string, bool>::iterator cur;
	cur=parameterUsed.begin();
	curParameter=mapParameter.begin();
	for(;cur!=parameterUsed.end() && curParameter!=endParameter;++cur, ++curParameter){
		if(!cur->second){
			if(parameterList!="") parameterList+=", ";
			parameterList+=curParameter->first;
		}
	}
	return parameterList;
}
