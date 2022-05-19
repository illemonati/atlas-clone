/*
 * TAllelicDepthCounts.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TALLELICDEPTHCOUNTS_H_
#define TALLELICDEPTHCOUNTS_H_

#include <stdint.h>
#include <string>

#include "TGenome.h"
#include "TTask.h"

namespace GenomeTasks{

//------------------------------------------
// TAllelicDepthCounts
//------------------------------------------
class TAllelicDepthCounts{
private:
	uint32_t _size            = 0;
	std::vector<uint32_t> _counts;

	void _allocateStorage(uint32_t MaxAllelicDepth);
public:
	TAllelicDepthCounts() = default;
	TAllelicDepthCounts(uint32_t MaxAllelicDepth);

	void resize(uint32_t MaxAllelicDepth);
	void clear();
	void addSite(const GenotypeLikelihoods::TBaseCounts & alleleCounts);
	void addSiteZeroDepth();
	void write(const std::string &filename, bool printEmpty);
};

//------------------------------------------
// TAllelicDepth
//------------------------------------------
class TAllelicDepth:public TGenome_windows{
private:
	TAllelicDepthCounts _counts;
	bool _writeEmpty;

	void _handleWindow() override;
	void _handleAlignment() override {}

public:
	TAllelicDepth();
	void quantifyAlleleicDepth();
};

//-------------------------------------------
// TTask_allelicDepth
//-------------------------------------------
class TTask_allelicDepth:public coretools::TTask{
public:
	TTask_allelicDepth(){ _explanation = "Writing genotype likelihoods to a GLF file"; };

	void run(){
		TAllelicDepth allelicDepth;
		allelicDepth.quantifyAlleleicDepth();
	}
};


}; // end namespace

#endif /* TALLELICDEPTHCOUNTS_H_ */
