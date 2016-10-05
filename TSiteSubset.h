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
#include "bamtools/utils/bamtools_fasta.h"
#include "TLog.h"

//store sites 0-based!

class TSiteSubsetWindow{
public:
	bool hasData;
	long start, end;
	std::map< long, std::pair<char,char> > positions; //stores reference and alternative allele
	std::map< int, char > positions2; //stores reference and alternative allele

	TSiteSubsetWindow(long Start, long End){
		hasData = false;
		start = Start;
		end = End;
	};
	~TSiteSubsetWindow(){};

	void addPosition(long pos, char & ref, char & alt){
		positions.insert(std::pair<long, std::pair<char,char> >(pos, std::pair<char,char>(ref, alt)));
	};

	void print(){
		std::cout << "[" << start << ", " << end << "]:";
		for(std::map<long,std::pair<char,char> >::iterator it=positions.begin(); it!=positions.end(); ++it) std::cout << " " << it->first << "(" << it->second.first << "," << it->second.second << ")";
		std::cout << std::endl;
	};

	long size(){
		return positions.size();
	};
};

class TSiteSubsetChr{
public:
	std::string name;
	std::map<int, TSiteSubsetWindow*> windows;
	std::map<int, TSiteSubsetWindow*>::iterator windowIt;
	int windowSize;
	int chrNumberInFasta;


	TSiteSubsetChr(std::string & Name, int & WindowSize){
		name = Name;
		windowSize = WindowSize;
		chrNumberInFasta = -1;
	};

	TSiteSubsetChr(std::string & Name, int & WindowSize, BamTools::SamHeader bamHeader){
		name = Name;
		windowSize = WindowSize;
		//find number by parsing through bam header
		int i=0;
		chrNumberInFasta = -1;
		for(BamTools::SamSequenceIterator chrIterator = bamHeader.Sequences.Begin(); chrIterator != bamHeader.Sequences.End(); ++chrIterator, ++i){
			if(chrIterator->Name == name) chrNumberInFasta = i;
		}
		if(chrNumberInFasta == -1) throw "chromosome '" + name + "' not present in FASTA refrence!";
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
	};

	void addPosition(std::vector<std::string> & tmp, const std::string & chr){
		long pos = stringToLong(tmp[1]) - 1; //make 0-based
		char ref = tmp[2][0];
		char alt = tmp[3][0];

		//check
		if(ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T'){
			std::string error = "Unknown reference allele '";
			error += ref;
			error += "' on chr " + chr;
			throw error + " at " + toString(pos) + "!";
		}
		if(alt != 'A' && alt != 'C' && alt != 'G' && alt != 'T'){
			std::string error = "Unknown alternative allele '";
			error += alt;
			error += "' on chr " + chr;
			throw error + " at " + toString(pos) + "!";
		}
		if(ref == alt) throw "Reference allele = alternative allele on chr " + chr + " at " + toString(pos+1) + "!";

		//identify window
		findOrCreateWindow(pos);
		windowIt->second->addPosition(pos, ref, alt);
	};

	bool addPosition(std::vector<std::string> & tmp, const std::string & chr, BamTools::Fasta & reference, std::string & error){
		long pos = stringToLong(tmp[1]) - 1; //make 0-based
		char ref = tmp[2][0];
		char alt = tmp[3][0];

		//check with reference
		char inRef;
		reference.GetBase(chrNumberInFasta, pos, inRef);
		if(ref != inRef){
			if(alt == inRef){
				//swap
				char tmp = alt;
				alt = ref;
				ref = tmp;
			} else {
				//problematic position -> skip!
				//logfile->warning("Conflict with FASTA reference, skipping position " + toString(pos+1) + " on chr " + name + "!");
				error = chr + "\t" + tmp[1] + "\t" + inRef + "\t" + ref + "\t" + alt;
				return false;
			}
		}

		//check
		if(ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T'){
			error = chr + "\t" + tmp[1] + "\t" + inRef + "\t" + ref + "\t" + alt;
			return false;
		}
		if(alt != 'A' && alt != 'C' && alt != 'G' && alt != 'T'){
			error = chr + "\t" + tmp[1] + "\t" + inRef + "\t" + ref + "\t" + alt;
			return false;
		}
		if(ref == alt){
			error = chr + "\t" + tmp[1] + "\t" + inRef + "\t" + ref + "\t" + alt;
			return false;
			//throw "Reference allele = alternative allele on chr " + chr + " at " + toString(pos+1) + "!";
		}

		//identify window
		findOrCreateWindow(pos);
		windowIt->second->addPosition(pos, ref, alt);
		return true;
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

	long size(){
		long s = 0;
		for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt)
			s += windowIt->second->size();
		return s;
	};
};

class TSiteSubset{
private:
	std::map<std::string, TSiteSubsetChr*> chromosomes;
	std::map<std::string, TSiteSubsetChr*>::iterator chrIt;
	int windowSize;
	std::string curChr;

	void readFile(TLog* logfile){
		logfile->listFlush("Reading sites to be used from '" + filename + "' ...");
		//open file
		std::ifstream sitesFile(filename.c_str());
		if(!sitesFile) throw "Failed to open sites file '" + filename + "'!";

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
				chrIt->second->addPosition(vec, chrIt->first);
			}
		}

		//close file
		sitesFile.close();

		//report
		long size = 0;
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
			size += chrIt->second->size();
		}
		logfile->write(" done!");
		logfile->conclude("Parsed " + toString(size) + " sites on " + toString(chromosomes.size()) + " chromosomes.");
	};

	void readFile(BamTools::Fasta & reference, BamTools::SamHeader bamHeader, TLog* logfile){ //version that checks witth fasta reference
		logfile->listFlush("Reading sites to be used from '" + filename + "' ...");
		//open file
		std::ifstream sitesFile(filename.c_str());
		if(!sitesFile) throw "Failed to open sites file '" + filename + "'!";

		//tmp variables
		long lineNum = 0;
		std::vector<std::string> vec;
		curChr = "";
		std::vector<std::string> conflictsWithReference;
		std::string error;

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
						chromosomes.insert(std::pair<std::string, TSiteSubsetChr*>(vec[0], new TSiteSubsetChr(vec[0], windowSize, bamHeader)));
						chrIt = chromosomes.find(vec[0]);
					}
					curChr = vec[0];
				}

				//add positions
				if(!chrIt->second->addPosition(vec, chrIt->first, reference, error)){
					//conflict with fasta -> add to vector
					conflictsWithReference.push_back(error);
				}
			}
		}
		//report
		logfile->write(" done!");
		long size = 0;
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
			size += chrIt->second->size();
		}
		logfile->conclude("Parsed " + toString(size) + " sites on " + toString(chromosomes.size()) + " chromosomes.");

		//write conflicts, if any
		if(conflictsWithReference.size() > 0){
			logfile->conclude("Reference conflicted with provided alleles at " + toString(conflictsWithReference.size()) + " positions!");

			std::string conflictFileName = filename + ".conflicts";
			std::ofstream out(conflictFileName.c_str());
			out << "chr\tpos\tref\tallele1\tallele2\n";
			for(std::vector<std::string>::iterator it = conflictsWithReference.begin(); it != conflictsWithReference.end(); ++it)
				out << *it << "\n";
			out.close();
			logfile->conclude("These positions were written to " + conflictFileName);
		}

		//close file
		sitesFile.close();

	};

public:
	std::string filename;
	TSiteSubset(std::string Filename, int & WindowSize, TLog* logfile){
		filename = Filename;
		windowSize = WindowSize;
		readFile(logfile);
		curChr = "";
	};

	TSiteSubset(std::string Filename, BamTools::Fasta & reference, BamTools::SamHeader bamHeader, int & WindowSize, TLog* logfile){
		filename = Filename;
		windowSize = WindowSize;
		readFile(reference, bamHeader, logfile);
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
		if(chrIt == chromosomes.end()){
			return false;
		}
		else return chrIt->second->hasPositionsInWindow(windowStart);
	}

	std::map<long,std::pair<char,char> >& getPositionInWindow(long & windowStart){
		//find chromosome
		chrIt = chromosomes.find(curChr);
		if(chrIt == chromosomes.end()) throw "TSiteSubset Error: chromosome '" + curChr + "' does not exist!";
		return chrIt->second->getPositionInWindow(windowStart);
	};

	long size(){
		long size = 0;
		for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
			size += chrIt->second->size();
		return size;
	};

};


#endif /* TSITESUBSET_H_ */
