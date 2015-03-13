/*
 * loci.h
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#ifndef LOCI_H_
#define LOCI_H_

#include <math.h>
#include "stringFunctions.h"
#include "TLog.h"
#include "bamtools/api/BamReader.h"
#include "bamtools/api/SamSequenceDictionary.h"



class TSite{
public:
	double errorRate;
	int posThreePrime, posFivePrime;
	double* emissionProbabilities;

	TSite(double ErrorRate, int PosThreePrime, int PosFivePrime){
		errorRate = ErrorRate;
		posThreePrime = PosThreePrime;
		posFivePrime = PosFivePrime;

		//fill emission probabilities
		emissionProbabilities = new double[10];
		fillEmissionProbabilities();
	};

	virtual ~TSite(){
		delete[] emissionProbabilities;
	};

	void update(double ErrorRate, int PosThreePrime, int PosFivePrime){
		errorRate = ErrorRate;
		posThreePrime = PosThreePrime;
		posFivePrime = PosFivePrime;
	};

	virtual void fillEmissionProbabilities(){
		fillEmissionProbabilities();
	};

	virtual double getEmissionProbability(int genotype){
		return emissionProbabilities[genotype];
	}
};

class TLocus{
public:
	double pHet, pHom;
	long pos;

	TLocus(std::vector<std::string> & line){
		pos = stringToLong(line[1]);
		pHom = stringToDouble(line[2]);
		pHet = stringToDouble(line[3]);
	};

	void addToSumsForEM(double & sumHom, double & sumHet, double & Het, double & oneMinusHet){
		double tmp = pHom * oneMinusHet + pHet * Het;
		sumHom += pHom / tmp;
		sumHet += pHet / tmp;
	};

	double LogLikelihood(double & Het, double & oneMinusHet){
		return log(pHom * oneMinusHet + pHet * Het);
	};
};

class TGenome{
public:
	//std::vector<TChromosome> chromosomes;
	std::string filename;
	TLog* logfile;
	BamTools::BamReader bamReader;
 	BamTools::BamRegion bamRegion;
 	BamTools::SamHeader bamHeader;
 	BamTools::BamAlignment bamAlignement;

 	BamTools::SamSequenceIterator chrIterator;
 	int chrNumber;
 	long chrLength;
 	long curStart, curEnd;
 	long windowSize;

	TGenome(TLog* Logfile){
		filename = "";
		logfile = Logfile;
		curStart = -1;
		curEnd = -1;
		windowSize = -1;
		chrNumber = -1;
		chrLength = -1;
	};

	TGenome(TLog* Logfile, std::string Filename, long WindowSize){
		logfile = Logfile;
		filename = Filename;
		windowSize = WindowSize;
		chrNumber = -1;
		chrLength = -1;
		curStart = -1;
		curEnd = -1;

		//open BAM file
		if (!bamReader.Open(filename))
			throw "Failed to open BAM file '" + filename + "'!";

		//load index file
		if(!bamReader.LocateIndex())
			throw "No index file found for BAM file '" + filename + "'!";

		//read header
		bamHeader = bamReader.GetHeader();
		chrIterator = bamHeader.Sequences.End();



	};
/*
	bool nextChromosome(){
		if(chrIterator == bamHeader.Sequences.End()){
			chrIterator = bamHeader.Sequences.Begin();
			chrNumber = 0;
		}
		++chrIterator; ++chrNumber;


		chrLength = stringToLong(chrIterator->Length);
		curStart = 0;
		curEnd = windowSize;

		logfile->endIndent();
		logfile->startIndent("Parsing chromosome '" + chrIterator->Name + "':");

		return true;
	}
	*/

	bool nextChromosome(){
		if(chrIterator == bamHeader.Sequences.End()){
			chrIterator = bamHeader.Sequences.Begin();
			chrNumber = 0;
		} else logfile->endIndent();

		//move to next
		++chrIterator; ++chrNumber;

		//did we reach end?
		if(chrIterator == bamHeader.Sequences.End()){
			return false;
		}

		//restart windows
		chrLength = stringToLong(chrIterator->Length);
		curStart = 0;
		curEnd = 0;

		logfile->endIndent();
		logfile->startIndent("Parsing chromosome '" + chrIterator->Name + "':");

		return true;
	}

	bool nextWindow(){
		//move to next region
		curStart = curEnd;
		curEnd += windowSize;
		if(curEnd > chrLength) curEnd = chrLength + 1;

		if(curStart > chrLength) return false;
		return bamReader.SetRegion(chrNumber, curStart, chrNumber, curEnd);

		while(bamReader.GetNextAlignment(bamAlignement)){
			std::cout << bamAlignement.RefID << std::endl;
		}

	};

	bool readData(){
		logfile->listFlush("Reading data on '" + chrIterator->Name + "' at [" + toString(curStart) + ", " + toString(curEnd) + ") ...");

		//parse through reads
		int coverage = 0;
		coverage = bamReader.GetReferenceCount();
		/*
		while(bamReader.GetNextAlignment(bamAlignement)){



			if(!bamAlignement.IsFailedQC() && !bamAlignement.IsDuplicate()){
				++coverage;
			}
		}
*/
		logfile->write(" done!");
		logfile->conclude("coverage = " + toString(coverage));
		return true;
	};
/*
	void readProbabilities(const std::string Filename, long maxPos=-1){
		//clear current data
		chromosomes.clear();

		//open file
		filename = Filename;
		if(maxPos > 0) logfile->startIndent("Reading probabilities for up to " + toString(maxPos) + " loci from file '" + filename + "':");
		else logfile->startIndent("Reading probabilities from file '" + filename + "':");
		std::ifstream file(filename.c_str());
		if(!file) throw "Could not open file '" + filename + "!";

		//parse file
		std::vector<std::string> line;
		std::string curChr = "";
		std::vector<TChromosome>::reverse_iterator revChrIt;
		long counter = 0;

		while(file.good() && !file.eof()){
			fillVectorFromLineWhiteSpace(file, line);
			if(line.size() > 0){
				if(line.size() != 4) throw "Wrong number of columns (expect chr, pos, pHomo, pHet) on line:\n" + concatenateString(line, "\t");

				//check if we move to new chr
				if(line[0] != curChr){
					//report old
					if(curChr != ""){
						logfile->write(" done!");
						logfile->conclude("Read " + toString(revChrIt->size()) + " loci");
					}

					//move
					curChr = line[0];
					chromosomes.push_back(TChromosome(line[0]));
					revChrIt = chromosomes.rbegin();
					logfile->listFlush("Parsing data from chromosome '" + curChr + "' ...");
				}

				//append data to chromosome
				revChrIt->addLocus(line);

				//count
				++counter;
				if(counter == maxPos) break;
			}
		}

		//report last
		logfile->write(" done!");
		logfile->conclude("Read " + toString(revChrIt->size()) + " loci");
		logfile->endIndent();
	};

	double runEM(int numIter, double tol, double initHet){
		//run EM to find estimate of heterozygosity
		logfile->startIndent("Running EM algorithm to estimate heterozygosity:");

		//variables
		double newHet = 0.01;
		double curHet = 0.01;
		double oneMinusHet;
		double sumHom, sumHet;
		double LL = -999999999;
		double newLL = -999999999;
		std::string progressString;
		int chrCounter = 0;

		//run iterations
		for(int i=0; i<numIter; ++i){
			logfile->startIndent("Interation " + toString(i+1) + ":");
			logfile->list("Initial het = " + toString(curHet));
			progressString = "Updating het... ";
			logfile->listFlush(progressString + "(0 / " + toString(chromosomes.size()) + " chromosomes)");

			//Calculate required sums
			curHet = newHet;
			oneMinusHet = 1.0 - curHet;
			sumHom = 0.0; sumHet = 0.0;
			chrCounter = 0;

			for(chrIt = chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++chrCounter){
				chrIt->addToSumsForEM(sumHom, sumHet, curHet, oneMinusHet);
				logfile->listOverFlush(progressString + "(" + toString(chrCounter) + " / " + toString(chromosomes.size()) + " chromosomes)");
			}
			logfile->overList(progressString + "done!                    ");

			//update F
			newHet = curHet * sumHet / (curHet * sumHet + oneMinusHet * sumHom);
			logfile->conclude("new het = " + toString(newHet));
			oneMinusHet = 1.0 - newHet;

			//break?
			if(newHet == curHet){
				logfile->conclude("H did not change -> stopping EM.");
				break;
			}

			//calculate Likelihood
			progressString = "Calculating LL with het = " + toString(newHet) + " ... ";
			chrCounter = 0;
			LL = newLL;
			newLL = 0.0;
			for(chrIt = chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt, ++chrCounter){
				newLL += chrIt->LogLikelihood(newHet, oneMinusHet);
				logfile->listFlush(progressString + "(0 / " + toString(chromosomes.size()) + " chromosomes)");
			}
			logfile->overList(progressString + "done!                    ");
			logfile->conclude("LL = " + toString(newLL));
			if(i > 0){
				double improvement = newLL - LL;
				logfile->conclude("LL improvement = " + toString(improvement));
				if(improvement < tol){
					logfile->conclude("LL improvement " + toString(improvement) + " < " + toString(tol) + " -> stopping EM.");
					break;
				}
			}
			logfile->endIndent();
		}

		//conclude
		logfile->endIndent();
		logfile->startIndent("MLE estimate:");
		logfile->list("H = " + toString(newHet));
		logfile->list("LL = " + toString(newLL));
		logfile->endIndent();
		logfile->endIndent();
		return(newHet);
	};
*/

	void test(){
		while(nextChromosome()){
			while(nextWindow()){
				readData();
			}
		}
	}


};




#endif /* LOCI_H_ */
