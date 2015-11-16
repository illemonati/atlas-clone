/*
 * TSiteSubset.h
 *
 *  Created on: Nov 16, 2015
 *      Author: wegmannd
 */

#ifndef TSITESUBSET_H_
#define TSITESUBSET_H_

#include <fstream>
#include <vector>
#include <map>

class TSiteSubsetWindow{
public:
	bool hasData;
	long start, end;
	std::map<long,std::pair<char,char> > positions; //stores reference and alternative allele

	TSiteSubsetWindow(long Start, long End){
		hasData = false;
		start = Start;
		end = End;
	};
	~TSiteSubsetWindow(){};

	void addPosition(long & pos, char & ref, char & alt, std::string & chr){
		if(ref != 'A' || ref != 'C' || ref != 'G' || ref != 'T') throw "Unknown reference allele '" + ref + "' on chr " + chr + " at " + toString(pos) + "!";
		if(alt != 'A' || alt != 'C' || alt != 'G' || alt != 'T') throw "Unknown alternative allele '" + ref + "' on chr " + chr + " at " + toString(pos) + "!";
		if(ref == alt) throw "Reference allele = alternative allele on chr " + chr + " at " + toString(pos) + "!";
		positions.insert(pos, std::pair<char,char>(ref, alt));
	};

	void print(){
		std::cout << "[" << start << ", " << end << "]:";
		for(std::map<long,std::pair<char,char> >::iterator it=positions.begin(); it!=positions.end(); ++it) std::cout << " " << it->first << "(" << it->second.first << "," << it->second.second << ")";
		std::cout << std::endl;
	};
};

class TSiteSubsetChr{
public:
	std::string name;
	std::map<int, TSiteSubsetWindow*> windows;
	std::map<int, TSiteSubsetWindow*>::iterator windowIt;
	int windowSize;

	TSiteSubsetChr(std::string & Name, int & WindowSize){
		name = Name;
		windowSize = WindowSize;
	};

	~TSiteSubsetChr(){
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
			windows.insert(std::pair<int, TSiteSubsetWindow*>(w, new TSiteSubsetWindow(w*windowSize, (w+1)*windowSize - 1)));
			findWindow(pos);
		}
	}

	void addPosition(std::vector<std::string> & tmp){
		long pos = stringToLong(tmp[1]);
		char ref = tmp[2];
		char alt = tmp[3];

		//identify window
		findOrCreateWindow(pos);
		windowIt->second->addPosition(pos, ref, alt);
	};

	void print(){
		std::cout << "Chromosome '" << name << "':" << std::endl;
		for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt) windowIt->second->print();
	};

	bool hasPositionsInWindow(const long & windowStart){
		findWindow(windowStart);
		if(windowIt == windows.end()) return false;
		return true;
	};

	std::map<long,std::pair<char,char> >& getPositionInWindow(const long & windowStart){
		findWindow(windowStart);
		if(windowIt == windows.end()) throw "TSiteSubset Error: window '" + toString(windowStart) + "' does not exist!";
		return windowIt->second->positions;
	};

};

class TSiteSubset{
private:
	std::map<std::string, TSiteSubsetChr*> chromosomes;
	std::map<std::string, TSiteSubsetChr*>::iterator chrIt;
	int windowSize;
	std::string curChr;

	void readFile(){
		//open file
		std::ifstream sitesFile(filename.c_str());
		if(!sitesFile) throw "Failed to open BED file '" + filename + "'!";

		//tmp variables
		long lineNum = 0;
		std::vector<std::string> vec;
		curChr = "";

		//read file
		while(sitesFile.good() && !sitesFile.eof()){
			++lineNum;
			fillVectorFromLineWhiteSpaceSkipEmpty(sitesFile, vec);
			//skip empty lines
			if(vec.size() > 0){
				if(vec.size() != 4) throw "Wrong number of columns in sites file '" + filename + "' on line " + toString(lineNum) + "!";

				//get chromosome
				if(vec[0] != curChr){
					chrIt = chromosomes.find(vec[0]);
					if(chrIt == chromosomes.end()){
						chromosomes.insert(std::pair<std::string, TSiteSubsetChr*>(vec[0], new TSiteSubsetChr(vec[0], windowSize)));
						chrIt = chromosomes.find(vec[0]);
					}
					curChr = vec[0];
				}

				//add positions
				chrIt->second->addPosition(vec);
			}
		}

		//close file
		sitesFile.close();
	};

public:
	std::string filename;
	TSiteSubset(std::string Filename, int & WindowSize){
		filename = Filename;
		windowSize = WindowSize;
		readFile();
		curChr = "";
	};

	~TSiteSubset(){
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
		std::cout << "Sites File '" << filename << "':" << std::endl;
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
		if(chrIt == chromosomes.end()) throw "TSiteSubset Error: chromosome '" + curChr + "' does not exist!";
		return chrIt->second->getPositionInWindow(windowStart);
	};

};


#endif /* TSITESUBSET_H_ */
