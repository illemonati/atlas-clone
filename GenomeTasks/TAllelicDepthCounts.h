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
#include "coretools/Main/TTask.h"

namespace GenomeTasks{

//------------------------------------------
// TAllelicDepthCounts
//------------------------------------------
class TAllelicDepthCounts{
private:
	size_t _size            = 0;
	std::vector<size_t> _counts;

public:
	TAllelicDepthCounts() = default;
	TAllelicDepthCounts(size_t MaxAllelicDepth);

	void resize(size_t MaxAllelicDepth);
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
	void run();
};

}; // end namespace

#endif /* TALLELICDEPTHCOUNTS_H_ */
