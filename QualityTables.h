/*
 * QualityTables.h
 *
 *  Created on: Nov 30, 2017
 *      Author: phaentu
 */

#ifndef QUALITYTABLES_H_
#define QUALITYTABLES_H_

#include <fstream>

//---------------------------------------------------------------
//TQualityTable
//---------------------------------------------------------------
class TQualityTable{
private:
	int maxQ;
	int maxQPlusOne;
	long* counts;
	double* freqs;
	bool initialized;

	//tmp
	long sum;

public:
	TQualityTable(){
		initialized = false;
		maxQ = -1;
		maxQPlusOne = -1;
		counts = NULL;
		freqs = NULL;
		sum = 0;
	};

	TQualityTable(int MaxQ){
		init(MaxQ);
	};

	TQualityTable(TQualityTable && other):maxQ(0),maxQPlusOne(0),counts(nullptr),freqs(nullptr),initialized(false),sum(0){
		//copy from other
		maxQ = other.maxQ;
		maxQPlusOne = other.maxQPlusOne;
		counts = other.counts;
		freqs = other.freqs;
		initialized = other.initialized;
		sum = other.sum;

		//set other to default
		other.maxQ = 0;
		other.maxQPlusOne = 0;
		other.counts = nullptr;
		other.freqs = nullptr;
		other.initialized = false;
		other.sum = 0;
	};


	TQualityTable& operator=(TQualityTable && other){
		if(this != &other){
			//copy from other
			maxQ = other.maxQ;
			maxQPlusOne = other.maxQPlusOne;
			counts = other.counts;
			freqs = other.freqs;
			initialized = other.initialized;
			sum = other.sum;

			//set other to default
			other.maxQ = 0;
			other.maxQPlusOne = 0;
			other.counts = nullptr;
			other.freqs = nullptr;
			other.initialized = false;
			other.sum = 0;
		}

		return *this;
	};

	void init(int MaxQ){
		maxQ = MaxQ; //Note: quality is phred(error) + 33!
		maxQPlusOne = maxQ + 1;
		counts = new long[maxQPlusOne];
		for(int i=33; i<maxQPlusOne; ++i)
			counts[i] = 0;
		freqs = new double[maxQPlusOne];
		sum = 0;
		initialized = true;
	};

	~TQualityTable(){
		if(initialized){
			delete[] counts;
			delete[] freqs;
		}
	};

	int getMaxQ(){
		return maxQ;
	};

	void add(const int & qual){
		if(qual < maxQPlusOne)
			++counts[qual];
	};

	void add(int* qual, int & len){
		for(int i=0; i<len; ++i){
			if(qual[i] < maxQPlusOne)
				++counts[qual[i]];
		}
	};

	void add(TQualityTable & other){
		int otherMaxQ = other.getMaxQ();
		int m = std::min(maxQ, otherMaxQ) + 1;
		for(int i=33; i<m; ++i)
			counts[i] += other.at(i);
	};

	long at(int qual){
		return counts[qual];
	};

	void calcFrequencies(){
		sum = 0;

		for(int i=33; i<maxQPlusOne; ++i)
			sum += counts[i];

		for(int i=33; i<maxQPlusOne; ++i)
			freqs[i] = (double) counts[i] / (double) sum;
	};

	void write(const std::string & filename){
		//open output file
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";

		//write header
		out << "Quality\tQuality(char)\tCounts\tFrequencies\tCumulativeFrequencies\n";

		//calc frequencies
		calcFrequencies();

		double cumulFreq = 0.0;
		for(int i=33; i<maxQPlusOne; ++i){
			cumulFreq += freqs[i];
			out << i-33 << "\t" << (char) i << "\t" << counts[i] << "\t" << freqs[i] << "\t" << cumulFreq << "\n";
		}
	};
};

//---------------------------------------------------------------
//TQualityTransformTable
//---------------------------------------------------------------
class TQualityTransformTable{
public:
	int maxQInTable;
	int maxQInTablePlusOne;
	double** table; //old qual / new qual
	bool initialized;

	TQualityTransformTable(int maxPhredIntInTable){
		initialized = false;
		initialize(maxPhredIntInTable);
	};

	TQualityTransformTable(){
		initialized = false;
		maxQInTable = 0;
		maxQInTablePlusOne = 0;
		table = NULL;
	};

	~TQualityTransformTable(){
		for(int i=0; i<maxQInTablePlusOne; ++i){
			delete[] table[i];
		}
		delete[] table;
	};

	void initialize(int maxPhredIntInTable){
		if(initialized == true)
			throw "Quality table already initialized!";

		maxQInTable = maxPhredIntInTable + 33;
		maxQInTablePlusOne = maxQInTable + 1;
		table = new double*[maxQInTablePlusOne];
		for(int i=0; i<maxQInTablePlusOne; ++i){
			table[i] = new double[maxQInTablePlusOne];
			for(int j=0; j<maxQInTablePlusOne; ++j){
				table[i][j] = 0;
			}
		}
		initialized = true;
	};

	void add(const int oldQuality, const int newQuality){
		if(oldQuality < maxQInTable && newQuality < maxQInTable){
			table[oldQuality][newQuality] += 1.0;
		}
	};

	double size(){
		double size = 0;
		for(int i=33; i<maxQInTablePlusOne; ++i){
			for(int j=33; j<maxQInTablePlusOne; ++j){
				size += table[i][j];
			}
		}
		return size;
	};

	void printTable(const std::string filename){
		//open file
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open output file '" + filename + "'!";

		//print header
		out << "oldQ/newQ";
		for(int i=33; i<maxQInTablePlusOne; ++i)
			out << "\t" << i-33;
		out << "\n";

		//get total
		double sum = size();

		//print rows
		for(int i=33; i<maxQInTablePlusOne; ++i){
			out << i-33;
			for(int j=33; j<maxQInTablePlusOne; ++j){
				out << "\t" << table[i][j] / sum;
			}
			out << "\n";
		}

		//close file
		out.close();
	};
};

//---------------------------------------------------------------
//TQualityTransformTables
//---------------------------------------------------------------
class TQualityTransformTables{
private:
	TReadGroups* readGroups;
	TQualityTransformTable* perReadGroupTables;
	TQualityTransformTable combinedTable;

public:
	TQualityTransformTables(TReadGroups & ReadGroups, int MaxQ){
		readGroups = &ReadGroups;

		combinedTable.initialize(MaxQ);
		perReadGroupTables = new TQualityTransformTable[readGroups->size()];
		for(int i=0; i<readGroups->size(); i++)
			perReadGroupTables[i].initialize(MaxQ);
	};

	~TQualityTransformTables(){
		delete[] perReadGroupTables;
	};

	void add(const int readGroup, const int oldQuality, const int newQuality){
		perReadGroupTables[readGroup].add(oldQuality, newQuality);
		combinedTable.add(oldQuality, newQuality);
	};

	void writeTables(std::string outputName){
		//print tables for read groups
		for(int i=0; i<readGroups->size(); ++i)
			perReadGroupTables[i].printTable(outputName + "_" + readGroups->getName(i) + "_qualityTransformation.txt");

		combinedTable.printTable(outputName + "_total_qualityTransformation.txt");
	};
};

#endif /* QUALITYTABLES_H_ */
