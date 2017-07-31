/*
 * TGLF.h
 *
 *  Created on: Jul 23, 2017
 *      Author: phaentu
 */

#ifndef TGLF_H_
#define TGLF_H_

#include <zlib.h>
#include "gzstream.h"
#include <algorithm>
//#include <bitset>

class TGlfHandle{
protected:
	std::string filename;
	gzFile gzfp;
	bool isOpen;
	long positionInFile;
	uint32_t offset;
	std::string version;
	uint8_t zero8, one8;
	uint32_t zero32;
	std::string curChr;
	std::string header;

public:
	TGlfHandle(){
		isOpen = false;
		gzfp = NULL;
		positionInFile = 0;
		offset = 0;
		version = "GLF3";
		zero8 = 0;
		one8 = 1;
		zero32 = 0;
	};

	void close(){
		if(isOpen){
			gzclose(gzfp);
			isOpen = false;
		}
	};
};

class TGlfWriter:public TGlfHandle{
private:
	long oldPos;
	uint8_t recordType1;

	void writeHeader(){
		write(version.c_str(), 4*sizeof(char));

		if(header.length() > 0){
			uint32_t labelLength = header.size();
			write(&labelLength, sizeof(uint32_t));
			write(header.c_str(), labelLength * sizeof(char));
		} else
			write(&zero32, sizeof(uint32_t));
	};

	template <typename T>
	inline void write(const T* buf, size_t len){
		positionInFile += gzwrite(gzfp, buf, len);
	};

public:
	TGlfWriter(){
		init();
	};
	TGlfWriter(std::string Filename){
		init();
		open(Filename, "");
	};
	void init(){
		oldPos = 0;
		recordType1 = one8 << 4;
	};
	~TGlfWriter(){
		close();
	};

	//open & close streams
	void open(std::string Filename, std::string Header){
		filename = Filename;
		gzfp = NULL;
		gzfp = gzopen(filename.c_str(),"wb");

		if(gzfp == NULL)
			throw "Failed to open file '" + filename + "' for writing!";
		isOpen = true;
		positionInFile = 0;
		curChr = "";

		//write header
		header = Header;
		writeHeader();
	};

	void newChromosome(std::string name, uint32_t length){
		if(curChr != "")
			write(&zero8, sizeof(uint8_t));

		//write new chromosome: length of label, label, length of ref sequence
		uint32_t labelLength = name.size();

		write(&labelLength, sizeof(uint32_t));
		write(name.c_str(), name.length() * sizeof(char));
		write(&length, sizeof(uint32_t));

		//set oldPos and curChr
		oldPos = 0;
		curChr = name;
	};

	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, uint8_t* genoQualities, uint32_t & maxLL){
		//record type
		//TODO: add reference
		write(&recordType1, sizeof(uint8_t));

		//offset
		offset = pos - oldPos;
		oldPos = pos;
		write(&offset, sizeof(uint32_t));

		//maxLL (capped at 255) and depth as one uint32_t
		if(maxLL > 255) maxLL = 255;
		uint32_t tmp = maxLL << 24;
		tmp = tmp | depth;
		write(&tmp, sizeof(uint32_t));

		//root mean square of mapping qualities
		write(&RMS_mappingQual, sizeof(uint8_t));

		//genotype likelihoods
		//TODO: test if using more accuracy on qualities matters
		write(genoQualities, 10*sizeof(uint8_t));
	};

};

class TGlfReader:public TGlfHandle{
private:
	bool reachedEndOfChr;
	uint32_t offset;
	uint32_t depth_mask;
	uint32_t tmpInt32;
	uint8_t tmpInt8;
	uint8_t tmpInt8_10[10];
	int _lenRead;
	int _i;

	void init(){
		isOpen = false;
		reachedEndOfChr = true;
		offset = 0;
		position = 0;
		maxLL = 0;
		depth = 0;
		RMS_mappingQual = 0;
		chr = "";
		chrLength = 0;
		depth_mask = 0xFFFFFF;
		tmpInt32 = 0;
		tmpInt8 = 0;
		_lenRead = 0;
		_i = 0;
	};

	template <typename T>
	inline bool read(T* buf, size_t len){
		_lenRead = gzread(gzfp, buf, len);
		if(!_lenRead) return false;
		positionInFile += _lenRead;
		return true;
	};

	bool readChr(){
		//read chromosome info
		uint32_t len;
		if(!read(&len, sizeof(uint32_t)))
				return false;
		char* tmp = new char[len];
		read(tmp, len*sizeof(char));
		chr.assign(tmp, len);

		read(&tmpInt32, sizeof(uint32_t));
		chrLength = tmpInt32;
		position = 0;

		return true;
	};

public:
	int recordType;
	long position;
	std::string chr;
	long chrLength;
	int maxLL;
	int depth;
	int RMS_mappingQual;
	int genotypeQualities[10];

	TGlfReader(){
		init();
	};
	TGlfReader(std::string Filename){
		init();
		open(Filename);
	};
	~TGlfReader(){
		close();
	};

	//open file and parse header
	void open(std::string Filename){
		filename = Filename;
		gzfp = NULL;
		gzfp = gzopen(filename.c_str(),"rb");

		if(gzfp == NULL)
			throw "Failed to open file '" + filename + "' for reading!";
		isOpen = true;
		positionInFile = 0;
		curChr = "";

		//parse header
		//version
		char buffer[4];
		read(buffer, 4*sizeof(char));
		version = buffer;

		uint32_t len;
		read(&len, sizeof(uint32_t));

		if(len > 0){
			char* header = new char[len];
			read(&header, len*sizeof(char));
		}

		//read info of first chromosome
		readChr();
	};

	bool readNext(){
		//read record type
		if(!read(&tmpInt8, sizeof(uint8_t)))
			return false;
		recordType = tmpInt8 >> 4;

		if(recordType == 0){
			readChr();
			return readNext();
		} else if(recordType == 1){
			//read data of a single position
			//offset
			read(&offset, sizeof(uint32_t));
			position += offset;

			//maxLL and depth
			read(&tmpInt32, sizeof(uint32_t));
			depth = tmpInt32 & depth_mask;
			maxLL = tmpInt32 >> 24;

			//root mean square of mapping qualities
			read(&tmpInt8, sizeof(uint8_t));
			RMS_mappingQual = (int) tmpInt8;

			//genotype likelihoods
			read(tmpInt8_10, 10*sizeof(uint8_t));

			for(_i=0; _i<10; ++_i){
				genotypeQualities[_i] = (int) tmpInt8_10[_i];
			}

			return true;
		} else throw "Unknown record type in file '" + filename + "'!";
	};

	void printChr(){
		std::cout << "CHROMOSOME: '" << chr << "' of length " << chrLength << "\n";
	};

	void printSite(){
		std::cout << chr << "\t" << position << "\t" << maxLL << "\t" << depth << "\t" << RMS_mappingQual;
		for(int i=0; i<10; ++i)
			std::cout << "\t" << genotypeQualities[i];
		std::cout << "\n";
	};

	void printToEnd(){ //For debugging
		//first print header
		std::cout << version << "\n";
		std::cout << header << "\n";
		printChr();

		//now parse file
		std::string oldChr = chr;
		while(readNext()){
			if(oldChr != chr){
				printChr();
				oldChr = chr;
			}
			printSite();
		}
	};

};

#endif /* TGLF_H_ */
