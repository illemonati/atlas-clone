/*
 * TAllelicDepthCounts.h
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#ifndef TALLELICDEPTHCOUNTS_H_
#define TALLELICDEPTHCOUNTS_H_

#include "TFile.h"
#include <map>

class TAllelicDepthCounts{
private:
	uint32_t _maxAllelicDepth;
	uint32_t _size;
	uint32_t**** _counts;


	void _allocateStorage(const uint32_t MaxAllelicDepth);
	void _freeStorage();

public:

	TAllelicDepthCounts(const uint32_t MaxAllelicDepth){
		_allocateStorage(MaxAllelicDepth);
	};
	~TAllelicDepthCounts(){
		_freeStorage();
	};

	void clear();
	void addSite(uint32_t numA, uint32_t numC, uint32_t numG, uint32_t numT);
	void addSiteZeroDepth();
	void write(const std::string filename, bool printEmpty);
};






#endif /* TALLELICDEPTHCOUNTS_H_ */
