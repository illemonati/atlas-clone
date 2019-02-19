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
	int i;
	long sum;

public:
	TQualityTable(){
		initialized = false;
		maxQ = -1;
		maxQPlusOne = -1;
		counts = NULL;
		freqs = NULL;
		i = 0;
		sum = 0;
	};

	TQualityTable(int MaxQ){
		init(MaxQ);
	};

	void init(int MaxQ){
		maxQ = MaxQ; //Note: quality is phred(error) + 33!
		maxQPlusOne = maxQ + 1;
		counts = new long[maxQPlusOne];
		for(i=33; i<maxQPlusOne; ++i)
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
		for(i=0; i<len; ++i){
			if(qual[i] < maxQPlusOne)
				++counts[qual[i]];
		}
	};

	void add(TQualityTable & other){
		int otherMaxQ = other.getMaxQ();
		int m = std::min(maxQ, otherMaxQ) + 1;
		for(i=33; i<m; ++i)
			counts[i] += other.at(i);
	};

	long at(int qual){
		return counts[qual];
	};

	void calcFrequencies(){
		sum = 0;

		for(i=33; i<maxQPlusOne; ++i)
			sum += counts[i];

		for(i=33; i<maxQPlusOne; ++i)
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
		for(i=33; i<maxQPlusOne; ++i){
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
	int maxQ;
	int maxQPlusOne;
	double** table; //old qual / new qual

	TQualityTransformTable(int MaxQ){
		maxQ = MaxQ + 33;
		maxQPlusOne = maxQ + 1;
		table = new double*[maxQPlusOne];
		for(int i=0; i<maxQPlusOne; ++i){
			table[i] = new double[maxQPlusOne];
			for(int j=0; j<maxQPlusOne; ++j){
				table[i][j] = 0;
			}
		}
	};

	~TQualityTransformTable(){
		for(int i=0; i<maxQPlusOne; ++i){
			delete[] table[i];
		}
		delete[] table;
	};

	void add(int oldQuality, int newQuality){
		if(oldQuality < maxQ && newQuality < maxQ){
			table[oldQuality][newQuality] += 1.0;
		}
	};

	double size(){
		double size = 0;
		for(int i=33; i<maxQPlusOne; ++i){
			for(int j=33; j<maxQPlusOne; ++j){
				size += table[i][j];
			}
		}
		return size;
	};

	void printTable(std::ofstream & out){
		//print header
		out << "oldQ/newQ";
		for(int i=33; i<maxQPlusOne; ++i)
			out << "\t" << i-33;
		out << "\n";

		//get total
		double sum = size();

		//print rows
		for(int i=33; i<maxQPlusOne; ++i){
			out << i-33;
			for(int j=33; j<maxQPlusOne; ++j){
				out << "\t" << table[i][j] / sum;
			}
			out << "\n";
		}
	};
};



#endif /* QUALITYTABLES_H_ */
