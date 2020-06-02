/*
 * TSoftClipping.h
 *
 *  Created on: Dec 9, 2019
 *      Author: wegmannd
 */

#ifndef TSOFTCLIPPING_H_
#define TSOFTCLIPPING_H_

#include "TGenome.h"

//TODO: does using "left" and "right" make sense? Should we not rather use 3' and 5' ends?

//--------------------------------------------------------
// TSoftClippingStatsFile
//--------------------------------------------------------
class TSoftClippingStatsFile{
private:
	TOutputFile _out;
	bool _printSequences;

public:
	TSoftClippingStatsFile(){ _printSequences = false; };
	TSoftClippingStatsFile(const std::string filename, const bool PrintSequences);
	void open(const std::string filename, const bool PrintSequences);
	void write(const BAM::TBamFile & bamFile);
};

//--------------------------------------------------------
// TAssessSoftClipping
//--------------------------------------------------------
class TAssessSoftClipping:public TGenome_filtered{
private:
	bool _writeAlignments;
	bool _printAll;

	TSoftClippingStatsFile statFile;

public:
	TAssessSoftClipping(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void assess();
};

//--------------------------------------------------------
// TRemoveSoftClippedBases
//--------------------------------------------------------
class TRemoveSoftClippedBases:public TGenome_filtered{
public:
	TRemoveSoftClippedBases(TParameters & Params, TLog* Logfile, TRandomGenerator* RandomGenerator);
	void removeSoftclippedBases();
};



#endif /* TSOFTCLIPPING_H_ */
