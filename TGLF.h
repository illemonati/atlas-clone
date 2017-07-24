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

class TGlfHandle{
protected:
	std::string filename;
	gzFile gzfp;
	bool isOpen;
	long positionInFile;
	std::string version;
	uint8_t zero, one;

public:
	TGlfHandle(){
		isOpen = false;
		gzfp = NULL;
		positionInFile = 0;
		version = "GLF3";
		zero = 0;
		one = 1;
	};
};

class TGlfWriter:public TGlfHandle{
private:
	long oldPos;
	std::string curChr;

public:
	TGlfWriter(){
		init();
	};
	TGlfWriter(std::string Filename){
		init();
		open(Filename);
	};
	void init(){
		oldPos = 0;
	};
	~TGlfWriter(){
		close();
	};

	//open & close streams
	void open(std::string Filename){
		filename = Filename;
		gzfp = NULL;
		gzFile gzfp = gzopen("out.gz","wb");
		if(gzfp == NULL)
			throw "Failed to open file '" + filename + "' for reading!";
		isOpen = true;
		positionInFile = 0;
		curChr = "";

		//write header
		writeHeader();
	};

	void writeHeader(){
		write(version.c_str(), 4*sizeof(char));
		write(0, sizeof(int32_t));
	};

	inline void write(const void* buf, size_t len){
		positionInFile += gzwrite(gzfp, buf, len);
	};

	void close(){
		if(isOpen)
			gzclose(gzfp);
	};

	void newChromosome(std::string name, long length){
		if(curChr != "")
			write(&zero, sizeof(uint8_t));

		//write new chromosome: length of label, label, length of ref sequence
		int32_t labelLength = name.size();
		write(&labelLength, sizeof(int32_t));
		write(name.c_str(), name.length() * sizeof(char));
		write(&length, sizeof(uint32_t));

		//set oldPos
		oldPos = 0;
	};

	void writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual, uint8_t* genoQualities){
		//record type
		write(&one, sizeof(uint8_t));

		//depth, capped at 255
		if(depth > 255) depth = 255;
		write(&depth, sizeof(uint32_t));

		//root mean square of mapping qualities
		write(&RMS_mappingQual, sizeof(uint8_t));

		//genotype likelihoods
		//TODO: test if using more accuracy on qualities matters
		write(genoQualities, 10*sizeof(uint8_t));
	};

};

/*
class TGlfReader:public TGlfHandle{
private:
	gz::igzstream handle;

public:
	TGlfReader(){
		isOpen = false;
	};
	TGlfReader(std::sting Filename){
		open(Filename);
	};
	~TGlfReader(){
		close();
	};

	//open & close streams
	void open(std::string Filename){
		filename = Filename;
		handle.open(filename.c_str());
		if(!handle)
			throw "Failed to open file '" + filename + "' for writing!";
		isOpen = true;

		//read header
	};

	void close(){
		if(isOpen)
			handle.close();
	};

};
*/


#endif /* TGLF_H_ */
