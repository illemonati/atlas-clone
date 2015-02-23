//---------------------------------------------------------------------------

#ifndef TParametersH
#define TParametersH
#include <map>
#include <vector>
#include <fstream>
#include "TLog.h"
#include "stringFunctions.h"
//---------------------------------------------------------------------------

class TParameters{
	private:
		std::string getParameter(std::string & my_name, bool mandatory);
		TLog* logfile;
		bool inputFileRead;
		std::string inputFileName;

	public:
		std::map<std::string, std::string> mapParameter;
		std::map<std::string, bool> parameterUsed;
		std::map<std::string, std::string>::iterator curParameter, endParameter;

		~TParameters(){ mapParameter.clear(); parameterUsed.clear(); };
		TParameters(std::vector<std::string> & commandLineParams, TLog* Logfile);
		TParameters(int & argc, char** argv, TLog* Logfile);
		void initialize(std::vector<std::string> & commandLineParams);
		void readInputfile(std::string fileName);
		bool parameterExists(std::string my_name);
		std::string getParameterString(std::string my_name, bool mandatory=true);
		std::string getParameterStringWithDefault(std::string my_name, std::string def);
		double getParameterDouble(std::string my_name, bool mandatory=true);
		double getParameterDoubleWithDefault(std::string my_name, double def);
		int getParameterInt(std::string my_name, bool mandatory=true);
		int getParameterIntWithDefault(std::string my_name, int def);
		long getParameterLong(std::string my_name, bool mandatory=true);
		long getParameterLongWithDefault(std::string my_name, long def);
		std::string getInputFile(){ return inputFileName;};
		bool hasInputFileBeenRead(){return inputFileRead;};

		void fillParameterIntoVector(std::string my_name, std::vector<int> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<long> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<double> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<float> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<bool> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<std::string> & vec, char delim, bool mandatory=true);
		std::string getListOfUnusedParameters();
};
#endif
