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
#include <map>

//read sorted bed files window by window
//store all data in chr / window combinations using vectors
//Store all positions 0-based, as in TWindow

class TBedReaderWindow{
public:
	bool hasData;
	long start, end;
	std::vector<long> positions;

	TBedReaderWindow(long Start, long End){
		hasData = false;
		start = Start;
		end = End;
	};
	~TBedReaderWindow(){};
	void addPosition(long & pos){
		positions.push_back(pos);
	};

	void print(){
		std::cout << "[" << start+1 << ", " << end+1 << "]:";
		for(std::vector<long>::iterator it=positions.begin(); it!=positions.end(); ++it) std::cout << " " << *it + 1;
		std::cout << std::endl;
	};

	long size(){
		return positions.size();
	}
};

class TBedReaderChromosome{
public:
	std::string name;
	std::map<int, TBedReaderWindow*> windows;
	std::map<int, TBedReaderWindow*>::iterator windowIt;
	int windowSize;


	TBedReaderChromosome(std::string & Name, int & WindowSize){
		name = Name;
		windowSize = WindowSize;
	};

	~TBedReaderChromosome(){
		//delete all windows
		for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt){
			delete windowIt->second;
		}
		windows.clear();
	};

	void findWindow(const long & pos){
		int w = (double) pos / (double) windowSize;
		windowIt = windows.find(w);
	}

	void findOrCreateWindow(const long & pos){
		findWindow(pos);
		if(windowIt == windows.end()){
			//insert window
			int w = (double) pos / (double) windowSize;
			windows.insert(std::pair<int, TBedReaderWindow*>(w, new TBedReaderWindow(w*windowSize, (w+1)*windowSize - 1)));
			findWindow(pos);
		}
	}

	void addPosition(std::vector<std::string> & tmp){
		long start = stringToLong(tmp[1]);
		long end = stringToLong(tmp[2]);

		//identify window
		findOrCreateWindow(start);

		//add position to that window
		//Note BED is already 0 indexed
		for(long i=start; i<end; ++i){
			if(i >= windowIt->second->end) findOrCreateWindow(i);
			windowIt->second->addPosition(i);
		}
	};

	void print(){
		for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt) windowIt->second->print();
	};

	bool hasPositionsInWindow(const long & windowStart){
		findWindow(windowStart);
		if(windowIt == windows.end()) return false;
		return true;
	};

	std::vector<long>& getPositionInWindow(long windowStart){
		findWindow(windowStart);
		if(windowIt == windows.end()) throw "TBedReader Error: window '" + toString(windowStart) + "' does not exist!";
		return windowIt->second->positions;
	};

	long size(){
		long s = 0;
		for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt)
			s += windowIt->second->size();
		return s;
	};
};

class TBedReader{
private:
	std::map<std::string, TBedReaderChromosome*> chromosomes;
	std::map<std::string, TBedReaderChromosome*>::iterator chrIt;
	int windowSize;
	std::string curChr;

	void readFile(BamTools::SamSequenceDictionary & Sequences, 	TLog* logfile){
		//open file
		std::ifstream bedFile(filename.c_str());
		if(!bedFile) throw "Failed to open BED file '" + filename + "'!";

		//tmp variables
		long lineNum = 0;
		std::vector<std::string> vec;
		curChr = "";

		//read file
		while(bedFile.good() && !bedFile.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(bedFile, vec);
			//skip empty lines
			if(vec.size() > 0){
				if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";

				//get chromosome
				if(!Sequences.Contains(vec[0])) logfile->warning("Chromosome '" + vec[0] + "' from BED file is not present in the BAM header!");
				if(vec[0] != curChr){
					chrIt = chromosomes.find(vec[0]);
					if(chrIt == chromosomes.end()){
						chromosomes.insert(std::pair<std::string, TBedReaderChromosome*>(vec[0], new TBedReaderChromosome(vec[0], windowSize)));
						chrIt = chromosomes.find(vec[0]);
					}
					curChr = vec[0];
				}

				//add positions
				chrIt->second->addPosition(vec);
			}
		}

		//close file
		bedFile.close();
	};

public:
	std::string filename;

	TBedReader(std::string Filename, int & WindowSize, BamTools::SamSequenceDictionary & Sequences, TLog* logfile){
		filename = Filename;
		windowSize = WindowSize;
		readFile(Sequences, logfile);
		curChr = "";
	};

	~TBedReader(){
		//delete all chromosomes
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
			delete chrIt->second;
		}
		chromosomes.clear();
	};

	void setChr(const std::string & chr){
		curChr = chr;
	};

	void print(){
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt) chrIt->second->print();
	};

	bool hasPositionsInWindow(const long & windowStart){
		chrIt = chromosomes.find(curChr);
		if(chrIt == chromosomes.end()) return false;
		else return chrIt->second->hasPositionsInWindow(windowStart);
	}

	std::vector<long>& getPositionInWindow(long & windowStart){
		//find chromosome
		chrIt = chromosomes.find(curChr);
		if(chrIt == chromosomes.end()) throw "TBedReader Error: chromosome '" + curChr + "' does not exist!";
		return chrIt->second->getPositionInWindow(windowStart);
	};

	long size(){
		long s=0;
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
			s += chrIt->second->size();
		return s;
	};

	int getNumChromosomes(){
		return chromosomes.size();
	};

};


#endif /* TBEDREADER_H_ */
