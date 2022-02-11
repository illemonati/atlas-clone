/*
 * TAllelicDepthCounts.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TALLELICDEPTHCOUNTS_H_
#define TALLELICDEPTHCOUNTS_H_

#include "TGenome.h"
#include <map>
#include "TTask.h"

namespace GenomeTasks{

//------------------------------------------
// TAllelicDepthCounts
//------------------------------------------
class TAllelicDepthCounts{
private:
	uint32_t _maxAllelicDepth;
	uint32_t _size;
	uint32_t**** _counts;
	bool _initialized;


	void _allocateStorage(const uint32_t MaxAllelicDepth);
	void _freeStorage();

public:
	TAllelicDepthCounts();
	TAllelicDepthCounts(const uint32_t MaxAllelicDepth);
	~TAllelicDepthCounts();

	void resize(const uint32_t MaxAllelicDepth);
	void clear();
	void addSite(const GenotypeLikelihoods::TBaseCounts & alleleCounts);
	void addSiteZeroDepth();
	void write(const std::string filename, bool printEmpty);
};

//------------------------------------------
// TAllelicDepth
//------------------------------------------
class TAllelicDepth:public TGenome_windows{
private:
	TAllelicDepthCounts _counts;
	bool _writeEmpty;

	void _handleWindow();

public:
	TAllelicDepth(coretools::TParameters & Parameters, coretools::TLog* Logfile, coretools::TRandomGenerator* RandomGenerator);
	void quantifyAlleleicDepth();
};

//-------------------------------------------
// TTask_allelicDepth
//-------------------------------------------
class TTask_allelicDepth:public coretools::TTask{
public:
	TTask_allelicDepth(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
		TAllelicDepth allelicDepth(Parameters, Logfile, _randomGenerator);
		allelicDepth.quantifyAlleleicDepth();
	}
};


}; // end namespace

#endif /* TALLELICDEPTHCOUNTS_H_ */
