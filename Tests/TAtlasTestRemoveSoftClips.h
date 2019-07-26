/*
 * TAtlasTestRemoveSoftClips.h
 *
 *  Created on: Jul 26, 2019
 *      Author: linkv
 */

#ifndef TESTS_TATLASTESTREMOVESOFTCLIPS_H_
#define TESTS_TATLASTESTREMOVESOFTCLIPS_H_

#include "TAtlasTest.h"
#include "TTaskList.h"
#include "bamtools/api/BamWriter.h"
#include "bamtools/api/SamHeader.h"
#include "bamtools/api/BamAlignment.h"

//------------------------------------------
//TAtlasTest_RemoveSoftClips
//------------------------------------------

class TAtlasTest_removeSoftClips:public TAtlasTest{
private:
	std::string filenameTag;
	std::string bamFileName;

	void writeBAM();
	bool checkNewBAM();

public:
	TAtlasTest_removeSoftClips(TParameters & params, TLog* logfile);
	~TAtlasTest_removeSoftClips(){};
	bool run();
};



#endif /* TESTS_TATLASTESTREMOVESOFTCLIPS_H_ */
