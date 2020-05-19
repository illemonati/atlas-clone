/*
 * TContextStats.h
 *
 *  Created on: Jan 28, 2020
 *      Author: wegmannd
 */

#ifndef TCONTEXTSTATS_H_
#define TCONTEXTSTATS_H_

#include <TBase.h>
#include "TGenotypeMap.h"
#include "TQualityMap.h"

class TContextStats{
private:
	unsigned int maxPhredInt;
	TGenotypeMap genoMap;
	TQualityMap qualMap;
	unsigned int** contextCounts;

	void _allocate(size_t MaxPhredInt){
		maxPhredInt = MaxPhredInt;

		//prepare table
	    contextCounts = new unsigned int*[maxPhredInt];
	    for(size_t i = 0; i < maxPhredInt; ++i){
	    	contextCounts[i] =  new unsigned int[genoMap.numContexts];
	    	for(size_t j=0; j<genoMap.numContexts; ++j)
	    		contextCounts[i][j] = 0;
	    }
	};

	void _free(){
		for(size_t i = 0; i < maxPhredInt; ++i){
		    	delete[] contextCounts[i];
		    }
		delete [] contextCounts;
	};

public:
	TContextStats(size_t MaxPhredInt){
		_allocate(MaxPhredInt);
	};

	~TContextStats(){
		_free();
	};

	void add(const TBase & base){
		++contextCounts[base.recalibratedQualityAsPhredInt][base.context];
	};

	void writeToFile(const std::string outputFileName){
		//open file
		TOutputFileZipped out(outputFileName);

		//write header
		std::vector<std::string> header = {"quality"};
		for(size_t j=0; j<genoMap.numContexts; ++j){
			header.push_back(genoMap.getGenotypeString(j));
		}
		out.writeHeader(header);

		//write data
		for(size_t i=0; i<maxPhredInt; ++i){
			out << i;
			for(size_t j=0; j<genoMap.numContexts; ++j){
				out << contextCounts[i][j];
			}
			out << std::endl;
		}
	};
};




#endif /* TCONTEXTSTATS_H_ */
