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
		//TLog* logfile;
		bool inputFileRead;
		std::string inputFileName;

		std::map<std::string, std::string> mapParameter;
		std::map<std::string, bool> parameterUsed;
		std::map<std::string, std::string>::iterator curParameter;

	public:
		~TParameters(){ mapParameter.clear(); parameterUsed.clear(); };
		TParameters();
		TParameters(std::vector<std::string> & commandLineParams, TLog* Logfile);
		TParameters(int & argc, char** argv, TLog* Logfile);
		void clear();

		//read parameters
		void addParameter(std::string name, std::string value);
		void initialize(std::vector<std::string> & commandLineParams, TLog* logfile);
		void readInputfile(std::string fileName, TLog* logfile);
		std::string getInputFile(){ return inputFileName;};
		bool hasInputFileBeenRead(){return inputFileRead;};

		//get single parameters
		bool parameterExists(std::string my_name);
		std::string getParameterString(std::string my_name, bool mandatory=true);
		double getParameterDouble(std::string my_name, bool mandatory=true);
		int getParameterInt(std::string my_name, bool mandatory=true);
		long getParameterLong(std::string my_name, bool mandatory=true);

		//get single parameters with Default
		std::string getParameterStringWithDefault(std::string my_name, std::string def);
		double getParameterDoubleWithDefault(std::string my_name, double def);
		int getParameterIntWithDefault(std::string my_name, int def);
		long getParameterLongWithDefault(std::string my_name, long def);

		//fill parameter vector
		void fillParameterIntoVector(std::string my_name, std::vector<int> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<long> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<double> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<float> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<bool> & vec, char delim, bool mandatory=true);
		void fillParameterIntoVector(std::string my_name, std::vector<std::string> & vec, char delim, bool mandatory=true);

		//fill parameter vector with default
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<int> & vec, char delim, std::string def);
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<long> & vec, char delim, std::string def);
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<double> & vec, char delim, std::string def);
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<float> & vec, char delim, std::string def);
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<bool> & vec, char delim, std::string def);
		void fillParameterIntoVectorWithDefault(std::string my_name, std::vector<std::string> & vec, char delim, std::string def);

		std::string getListOfUnusedParameters();
};
#endif
