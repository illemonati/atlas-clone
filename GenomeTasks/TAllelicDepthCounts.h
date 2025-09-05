/*
 * TAllelicDepthCounts.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TALLELICDEPTHCOUNTS_H_
#define TALLELICDEPTHCOUNTS_H_

#include <string>
#include <vector>

#include "TBamWindowTraverser.h"
#include "TSiteTraverser.h"
#include "genometools/Genotypes/Containers.h"

namespace GenomeTasks{

//------------------------------------------
// TAllelicDepthCounts
//------------------------------------------
class TAllelicDepthCounts{
private:
	size_t _size            = 0;
	std::vector<size_t> _counts;

public:
	void resize(size_t MaxAllelicDepth);
	void clear();
	void addSite(const genometools::TBaseCounts & alleleCounts);
	void write(const std::string &filename, bool printEmpty);
	size_t size() const noexcept {return _counts.size();}
};

//------------------------------------------
// TAllelicDepth
//------------------------------------------
class TAllelicDepth {
private:
	BAM::TSiteTraverser _siteTraverser;

	TAllelicDepthCounts _counts;
	bool _writeEmpty;

	void _traverseSites();

public:
	TAllelicDepth();
	void run();
};

}; // end namespace

#endif /* TALLELICDEPTHCOUNTS_H_ */
