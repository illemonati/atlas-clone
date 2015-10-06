/*
 * TBedReader.h
 *
 *  Created on: Oct 6, 2015
 *      Author: wegmannd
 */

#ifndef TBEDREADER_H_
#define TBEDREADER_H_

#include <fstream>
#include <vector>

//read sorted bed files window by window
//Note: it will never restart from the beginning

class TBedReader{
private:
	std::ifstream bedFile;
	bool fileOpen;
	std::vector<std::string> tmpVec;
	std::vector<long> positions;
	//storage for first base on next chromosome
	bool nextOnNewChr;
	bool storageUsed;
	std::string chr, nextChr;
	long nextPos, lastPos;

	void readLine(){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(bedFile, tmpVec);
		if(tmpVec.size() < 4) throw "Less than four columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";
		nextChr = tmpVec[0];
		lastPos = nextPos;
		nextPos = stringToLong(tmpVec[2]); //use third column because BED is 0 indexed!
		if(nextChr == chr && nextPos <= lastPos) throw "BED file '" + filename + "' seems not to be sorted!";
	};

public:
	std::string filename;
	long start, end;
	long lineNum;


	TBedReader(std::string Filename){
		filename = Filename;
		bedFile.open(filename.c_str());
		if(!bedFile) throw "Failed to open BED file '" + filename + "'!";
		fileOpen = true;
		start = 0;
		end = 0;
		lineNum = 0;
		chr = ""; nextChr = "";
		nextOnNewChr = true;
		storageUsed = false;
		nextPos = 0; lastPos = 0;
	};

	~TBedReader(){
		if(fileOpen) bedFile.close();
	};



	void read(long Size){
		//read up to size bases. Stop at end of chromosome.
		std::vector<std::string> vec;
		long pos;
		long readBases = 0;

		//check if we will be on a new chromosome
		if(nextOnNewChr){
			positions.clear();
			if(!storageUsed) readLine();
			else storageUsed = false;
			chr = nextChr;
			positions.push_back(nextPos);
			nextOnNewChr = false;
			++readBases;
		}

		//parse line by line
		while(readBases < Size && bedFile.good() && !bedFile.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(bedFile, vec);
			if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";

			//check if we jump to a new chromosome
			if(vec[0] != chr){
				storageUsed = true;
				break;
			}
			positions.push_back(nextPos);
		}

		//set start and end
		start = *positions.begin();
		end = *positions.rbegin();
	};


};


#endif /* TBEDREADER_H_ */
