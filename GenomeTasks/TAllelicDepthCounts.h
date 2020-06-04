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
	void addSite(uint32_t numA, uint32_t numC, uint32_t numG, uint32_t numT);
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
	TAllelicDepth(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void quantifyAlleleicDepth();
};




#endif /* TALLELICDEPTHCOUNTS_H_ */
