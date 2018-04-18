/*
 * TWindow.cpp
 *
 *  Created on: May 17, 2015
 *      Author: wegmannd
 */

#include "TWindow.h"


//-------------------------------------------------------
//Twindow
//-------------------------------------------------------
TWindow::TWindow(){
	start = -1;
	end = -1;
	length = -1;
	sites = NULL;
	sitesInitialized = false;
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesDepthAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
};

TWindow::TWindow(long Start, long End){
	start = Start;
	end = End;
	initSites(end - start); //end NOT in window!
};

void TWindow::initSites(long newLength){
	if(sitesInitialized){
		clear();
		delete[] sites;
	}
	length = newLength;
	if(length > 0){
		try{
			sites = new TSite[length];;
		} catch(...){
			throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
		}
	} else sites = NULL;

	sitesInitialized = true;
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesDepthAtLeastTwo = -1.0;
	numReadsInWindow = 0;
}

void TWindow::clear(){
	if(sitesInitialized){
		for(int i=0; i<length; ++i) sites[i].clear();
		delete[] sites;
	}
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesDepthAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
};

void TWindow::move(long Start, long End){
	start = Start;
	end = End;
	if(sitesInitialized){
		if((end - start) != length){
			initSites(end - start);
		} else {
			clear();
		}
	} else initSites(end - start);
};

bool TWindow::addFromRead(TAlignmentParser & alignmentParser, TPMD* pmdObjects){
	/* Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window
	 */

	//TODO: move to alignment parser instead?

	//check if alignment is inside window
	if(alignmentParser.position >= end) return true;
	if(alignmentParser.position + alignmentParser.length < start) return false;

	//find which position to consider first
	++numReadsInWindow;
	int firstPos = alignmentParser.position - start;

	//std::cout << "[" << start << "," << end <<  "]: firstPos = " << firstPos;

	int p = 0;

	if(firstPos < 0){
		while(p < alignmentParser.length && (firstPos + alignmentParser.alignedPos[p]) < 0)
			++p;
		if(p == alignmentParser.length)
			return false;
	}
	int internalPos;

	//std::cout << " -> firstPos = " << firstPos << ", p = " << p << std::endl;

	/* Note:
	 *  1) Reference is 5' -> 3'
	 *  2) distance is 0-based!
	 *  3) Ignoring indels in other mate when calculating distances
	 *  4) Function add needs first P(C->T), then P(G->A)
	 */

	for(; p < alignmentParser.length; ++p){
		if(alignmentParser.aligned[p] && alignmentParser.base[p] != N){
			internalPos = firstPos + alignmentParser.alignedPos[p];
			if(internalPos >= length)
				return true; //since part of the read maps to next window
			sites[internalPos].add(alignmentParser.base[p], alignmentParser.quality[p], p, alignmentParser.length-p, pmdObjects[alignmentParser.readGroupId].getProbCT(alignmentParser.distFrom5Prime[p]), pmdObjects[alignmentParser.readGroupId].getProbGA(alignmentParser.distFrom3Prime[p]), alignmentParser.context[p], alignmentParser.readGroupId);

			//std::cout << alignemntParser.position << "[" << p << "] -> " <<  << std::endl;
		}
	}

	return false;
}


void TWindow::addReferenceBaseToSites(BamTools::Fasta & reference, int & refId){
	if(!referenceBaseAdded){
		int stop = end - 1; //note that end is last position + 1
		std::string ref; //fasta object fills string
		reference.GetSequence(refId, start, stop, ref);
		for(int i=0; i<length; ++i){
	//		if(sites[i].hasData)
			sites[i].setRefBase(ref[i]);
		}
		referenceBaseAdded = true;
	}
}

void TWindow::addReferenceBaseToSites(TSiteSubset* subset){
	if(!referenceBaseAdded){
		if(subset->hasPositionsInWindow(start)){
			//now only run over sites listed in that window
			std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
			int pos;
			for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				pos = it->first - start;
				sites[pos].setRefBase(it->second.first);
			}
		}
		referenceBaseAdded = true;
	}
}

void TWindow::applyMask(TBedReader* mask, bool doInverseMasking){
	int pos;
	long first = start;
	if(doInverseMasking){
		if(mask->hasPositionsInWindow(start)){
			std::vector<long> thesePos = mask->getPositionInWindow(start);
			for(std::vector<long>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				pos = *it - start;
				//clear sites between regions (if there are none pos==first)
				for(int i=first; i<pos; ++i){
					if(pos < length){
						sites[i].clear();
					}
				}
				first = pos + 1;
			}
			//clear rest of window if necessary
			for(int i=first; i<end-start; ++i){
				sites[i].clear();
			}
		//else clear entire window
		} else	clear();
	} else{
		if(mask->hasPositionsInWindow(start)){
			std::vector<long> thesePos = mask->getPositionInWindow(start);
			//skip sites listed in mask by setting their hasData = false
			for(std::vector<long>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
				pos = *it - start;
				if(pos < length) sites[pos].clear();
			}
		}
	}
}

void TWindow::maskCpG(BamTools::Fasta & reference, int & refId){
	std::string ref; //fasta object fills string
	//note that end is last position + 1
	for(int i=0; i<length; ++i){
		if(ref[i+1] == 'C' && ref[i+2] == 'G') sites[i].clear();
		else if(ref[i] == 'C' && ref[i+1] == 'G') sites[i].clear();
	}
}

void TWindow::estimateBaseFrequencies(){
	//estimate initial base frequencies
	baseFreq.clear();
	for(int i=0; i<length; ++i){
		sites[i].addToBaseFrequencies(baseFreq);
	}
	baseFreq.normalize();
}

void TWindow::calculateEmissionProbabilities(TRecalibration* recalObject){
	for(int i=0; i<length; ++i){
		//std::cout << start + i << ": ";

		if(sites[i].hasData)
			recalObject->calcEmissionProbabilities(sites[i]);

		//std::cout << std::endl;
	}
}

void TWindow::callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef){
	if(isVCF){
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				if(sites[i].hasData) recalObject->calcEmissionProbabilities(sites[i]);
				std::string basesString = sites[i].getBases();
				sites[i].callMLEGenotypeVCF(genoMap, randomGenerator, out, gVCF, noAltIfHomoRef, basesString);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					recalObject->calcEmissionProbabilities(sites[i]);
					std::string basesString = sites[i].getBases();
					sites[i].callMLEGenotypeVCF(genoMap, randomGenerator, out, gVCF, noAltIfHomoRef, basesString);
					out << "\n";
				}
			}
		}
	} else {
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				if(sites[i].hasData) recalObject->calcEmissionProbabilities(sites[i]);
				sites[i].callMLEGenotype(genoMap, randomGenerator, out);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					recalObject->calcEmissionProbabilities(sites[i]);
					sites[i].callMLEGenotype(genoMap, randomGenerator, out);
					out << "\n";
				}
			}
		}
	}
}

void TWindow::printPileup(TRecalibration* recalObject, gz::ogzstream & out, std::string & chr){
	//print pileup
	for(int i=0; i<length; ++i){
		recalObject->calcEmissionProbabilities(sites[i]);
		out << chr << "\t" << start + i + 1;
		sites[i].printPileup(out);
		out << "\n";
	}
}

void TWindow::calcDepth(){
	//calculate and return coverage
	depth = 0.0;
	long noData = 0;
	long plentyData = 0;
	int cov;
	for(int i=0; i<length; ++i){
		cov = sites[i].depth();
		depth += cov;
		if(cov == 0) ++noData;
		else if(cov > 1) ++ plentyData;
	}

	depth = depth / (double) length;
	numSitesWithData = length - noData;
	fractionSitesNoData = (double) noData / (double) length;
	fractionsitesDepthAtLeastTwo = (double) plentyData / (double) length;
}

void TWindow::calcFracN(){
	double numN = 0.0;
	for(int i=0; i<length; ++i)	if(sites[i].referenceBase == 'N') ++numN;
	fractionRefIsN = numN / (double) length;
}

void TWindow::calcDepthPerSite(long* siteDepth, size_t maxDepth){
	//calculate and return coverage
	depth = 0.0;
	long noData = 0;
	long plentyData = 0;
	size_t cov;

	for(int i=0; i<length; ++i){
		cov = sites[i].depth();
		if(cov <= maxDepth)	siteDepth[cov] += 1;
		else siteDepth[maxDepth + 1] += 1; //else it should be in the "greater than" bin

		if(cov == 0) ++ noData;
		else if(cov > 1) ++ plentyData;
	}

	depth = depth / (double) length;
	fractionSitesNoData = (double) noData / (double) length;
	fractionsitesDepthAtLeastTwo = (double) plentyData / (double) length;
}

void TWindow::printDepthPerSite(gz::ogzstream & out, std::string & chr){
	//print depth for each site to file
	for(int i=0; i<length; ++i){
		out << chr << "\t" << start + i + 1 << "\t" << sites[i].depth() << "\n";
	}
}

void TWindow::countAlleles(long**** siteImbalance, const unsigned int & maxCov){
	//calculate and return imbalance
	for(int i=0; i<length; ++i){
		if(sites[i].depth() <= maxCov && sites[i].depth() > 0)
			sites[i].countAlleles(siteImbalance);
		else if(sites[i].depth()  == 0){
			++siteImbalance[0][0][0][0];
		}

	}
}

void TWindow::applyDepthFilter(int minDepth, size_t maxDepth){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepth || sites[i].bases.size() > maxDepth)
				sites[i].clear();
		}
	}
}

void TWindow::createDepthMask(size_t minDepthForMask, size_t maxDepthForMask, std::ofstream & outputMaskFile, std::string & chr){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepthForMask || sites[i].bases.size() > maxDepthForMask){
				outputMaskFile << chr << "\t" << start + i << "\t" << start + i + 1 << "\n";
			}
		}
	}
}

void TWindow::addSitesToBQSR(TRecalibrationBQSR & bqsr, TLog* logfile){
	logfile->listFlush("Adding sites to BQSR ...");
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			bqsr.addSite(sites[i]);
		}
	}
	logfile->done();
}

void TWindow::addSitesToBQSR(TRecalibrationBQSR & bqsr, TSiteSubset* subset, TLog* logfile){
	logfile->listFlush("Adding sites to BQSR ...");
	//now only run over sites listed in that window
	std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData){
			sites[pos].setRefBase(it->second.second);
			bqsr.addSite(sites[pos]);
		}
	}
	logfile->done();

}

void TWindow::addSitesToQualityTransformTable(TRecalibration* recalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile){
	logfile->listFlush("Adding sites to quality transformation tables ...");
	std::vector<TBase*>::iterator it;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			for(it = sites[i].bases.begin(); it != sites[i].bases.end(); ++it){
				QTtables.at((*it)->readGroup)->add((*it)->quality, recalObject->getQualityFromBase(**it));
				QTtables.at(QTtables.size() - 1)->add((*it)->quality, recalObject->getQualityFromBase(**it));
			}
		}
	}
	logfile->done();
}

void TWindow::addSitesToQualityTransformTable(TRecalibration* recalObject, TRecalibration* otherRecalObject, std::vector<TQualityTransformTable*> & QTtables, TLog* logfile){
	logfile->listFlush("Adding sites to quality transformation tables ...");
	std::vector<TBase*>::iterator it;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			for(it = sites[i].bases.begin(); it != sites[i].bases.end(); ++it){
				QTtables.at((*it)->readGroup)->add(recalObject->getQualityFromBase(**it), otherRecalObject->getQualityFromBase(**it));
				QTtables.at(QTtables.size() - 1)->add(recalObject->getQualityFromBase(**it), otherRecalObject->getQualityFromBase(**it));
			}
		}
	}
	logfile->done();
}

void TWindow::addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile){
	logfile->listFlush("Adding sites to PMD tables ...");
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			//pmdTables.add(sites[i]);
		}
	}
	logfile->done();
}

//-------------------------------------------------------
//TwindowDiploid
//-------------------------------------------------------
void TWindow::addSitesToThetaEstimator(TRecalibration* recalObject, TThetaEstimator & estimator){
	//first calculate emission probabilities
	calculateEmissionProbabilities(recalObject);

	//now add sites
	addSitesToThetaEstimator(estimator);
}

void TWindow::addSitesToThetaEstimator(TThetaEstimator & estimator){
	//assumes that emission probabilities were calculated
	for(int i=0; i<length; ++i){
		estimator.add(sites[i]);
	}
}

void TWindow::callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool & isVCF, bool & noAltIfHomoRef, bool & beagle, bool & printOnlyGL){
	//check if we need to process this window
	if(subset->hasPositionsInWindow(start)){
		//now only run over sites listed in that window
		std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
		int pos;
		for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			pos = it->first - start;
			if(beagle) {
				if(sites[pos].hasData) recalObject->calcEmissionProbabilities(sites[pos]);
				sites[pos].callMLEGenotypeKnownAllelesBeagle(genoMap, randomGenerator, out, it->second.second, chr, pos, start, printOnlyGL);
			} else {
				out << chr << "\t" << it->first + 1;
				if(sites[pos].hasData) recalObject->calcEmissionProbabilities(sites[pos]);
				if(isVCF){
					std::string basesString = sites[pos].getBases();
					sites[pos].callMLEGenotypeVCFKnownAlleles(genoMap, randomGenerator, out, it->second.second, noAltIfHomoRef, basesString);
				} else sites[pos].callMLEGenotypeKnownAlleles(genoMap, randomGenerator, out, it->second.second);
				out << "\n";
			}
		}
	}
}

void TWindow::callBayesianGenotype(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF){
	//calc prior probabilities on Genotypes
	double* pGenotype = new double[10];
	estimator.fillPGenotype(pGenotype);

	//now call genotypes. Note: emission probabilities have already been calculated when estimating theta!
	if(isVCF){
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				sites[i].callBayesianGenotypeVCF(pGenotype, genoMap, randomGenerator, out);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					sites[i].callBayesianGenotypeVCF(pGenotype, genoMap, randomGenerator, out);
					out << "\n";
				}
			}
		}
	} else {
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				sites[i].callBayesianGenotype(pGenotype, genoMap, randomGenerator, out);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					sites[i].callBayesianGenotype(pGenotype, genoMap, randomGenerator, out);
					out << "\n";
				}
			}
		}
	}

	//clean up
	delete[] pGenotype;
}

void TWindow::callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF){
	//check if we need to process this window
	if(subset->hasPositionsInWindow(start)){
		//calc prior probabilities on Genotypes
		double* pGenotype = new double[10];
		estimator.fillPGenotype(pGenotype);

		//now only run over sites listed in that window
		std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
		int pos;
		for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			pos = it->first - start;
			out << chr << "\t" << it->first + 1;
			if(isVCF)
				sites[pos].callBayesianGenotypeVCFKnownAlleles(pGenotype, genoMap, randomGenerator, out, it->second.second);
			else
				sites[pos].callBayesianGenotypeKnownAlleles(pGenotype, genoMap, randomGenerator, out, it->second.second);
			out << "\n";
		}

		//clean up
		delete[] pGenotype;
	}
}

void TWindow::callAllelePresence(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool noAltIfHomoRef){
	//calc prior probabilities on Genotypes
	double* pGenotype = new double[10];
	estimator.fillPGenotype(pGenotype);

	//now call allele presence. Note: emission probabilities have already been calculated when estimating theta!
	if(isVCF){
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				std::string baseString = sites[i].getBases();
				sites[i].callAllelePresenceVCF(pGenotype, genoMap, randomGenerator, out, noAltIfHomoRef, baseString);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					std::string baseString = sites[i].getBases();
					sites[i].callAllelePresenceVCF(pGenotype, genoMap, randomGenerator, out, noAltIfHomoRef, baseString);
					out << "\n";
				}
			}
		}
	} else {
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				sites[i].callAllelePresence(pGenotype, genoMap, randomGenerator, out);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					sites[i].callAllelePresence(pGenotype, genoMap, randomGenerator, out);
					out << "\n";
				}
			}
		}
	}

	//clean up
	delete[] pGenotype;
}

void TWindow::callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll){
	//now call allele presence. Note: emission probabilities have already been calculated when estimating theta!
	if(printAll){
		for(int i=0; i<length; ++i){
			out << chr << "\t" << start + i + 1;
			sites[i].callRandomBase(randomGenerator, out);
			out << "\n";
		}
	} else {
		for(int i=0; i<length; ++i){
			if(sites[i].hasData){
				out << chr << "\t" << start + i + 1;
				sites[i].callRandomBase(randomGenerator, out);
				out << "\n";
			}
		}
	}
}

void TWindow::majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll){
	//now call allele presence. Note: emission probabilities have already been calculated when estimating theta!
	if(printAll){
		for(int i=0; i<length; ++i){
			out << chr << "\t" << start + i + 1;
			sites[i].majorityCall(randomGenerator, out);
			out << "\n";
		}
	} else {
		for(int i=0; i<length; ++i){
			if(sites[i].hasData){
				out << chr << "\t" << start + i + 1;
				sites[i].majorityCall(randomGenerator, out);
				out << "\n";
			}
		}
	}
}

void TWindow::callAllelePresenceKnwonAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF, bool noAltIfHomoRef){
	//check if we need to process this window
	if(subset->hasPositionsInWindow(start)){
		//calc prior probabilities on Genotypes
		double* pGenotype = new double[10];
		estimator.fillPGenotype(pGenotype);

		//now only run over sites listed in that window
		std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
		int pos;
		for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			pos = it->first - start;
			out << chr << "\t" << it->first + 1;
			if(isVCF){
				std::string basesString = sites[pos].getBases();
				sites[pos].callAllelePresenceVCFKnownAlleles(pGenotype, genoMap, randomGenerator, out, it->second.second, noAltIfHomoRef, basesString);
			} else {
				sites[pos].callAllelePresenceKnownAlleles(pGenotype, genoMap, randomGenerator, out, it->second.second);
			}
			out << "\n";

		}

		//clean up
		delete[] pGenotype;
	}
}

void TWindow::addToGLF(TGlfWriter & writer, bool printAll){
	//TODO: calculate root mean squared mapping qualities for sites (now just passing 0). Would be helpful in VCFs as well
	uint8_t* gl = new uint8_t[10];
	uint32_t maxLL;
	if(printAll){
		for(int i=0; i<length; ++i){
			sites[i].calculateNormalizedGenotypeLikelihoods(gl, maxLL);
			writer.writeSite(start + i + 1, sites[i].depth(), 0, gl, maxLL);
		}
	} else {
		for(int i=0; i<length; ++i){
			if(sites[i].hasData){
				sites[i].calculateNormalizedGenotypeLikelihoods(gl, maxLL);
				writer.writeSite(start + i + 1, sites[i].depth(), 0, gl, maxLL);
			}
		}
	}
	delete[] gl;
}

void TWindow::generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine){
	//calc prior probabilities on Genotypes
	double* pGenotype = new double[10];
	estimator.fillPGenotype(pGenotype);

	//now call heterozygosity in blocks
	int nBlocks = length / blockSize;
	int start;
	double logPHomo;
	double logConfidence = log(confidence);
	double logConfidenceHet = log(1.0 - confidence);
	double tmp;

	//loop over blocks
	for(int b=0; b<nBlocks; ++b){
		start = b*blockSize;
		logPHomo = 0.0;

		for(int i=0; i<blockSize; ++i){
			if(sites[start + i].hasData){
				tmp = sites[start + i].calculatePHomozygous(pGenotype);
				logPHomo += log(tmp);
			}
		}

		//check if we are heterozygous
		if(logPHomo > logConfidence){
			out << 'T';
		} else if(logPHomo < logConfidenceHet){
			out << 'K';
		} else {
			out << 'N';
		}

		//do we add a new line?
		if(nCharOnLine == 59){
			nCharOnLine = 0;
			out << '\n';
		} else ++nCharOnLine;
	}
	delete[] pGenotype;
}


void TWindow::fillPGenotype(double* pGenotype){
	for(int i=0; i<4; ++i){
		pGenotype[i] = baseFreq[i];
	}
}

double TWindow::calcLogLikelihood(){
	double pGenotype[4];
	fillPGenotype(pGenotype);

	double LL = 0.0;
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			LL += sites[i].calculateLogLikelihood(pGenotype);
		}
	}
	return LL;
}

void TWindow::addToRecalibrationEM(TRecalibrationEM & recalObject){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject.addSite(sites[i]);
		}
	}
}

void TWindow::addToRecalibrationEM(TRecalibrationEM & recalObject, TSiteSubset* subset){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	//now only run over sites listed in that window
	std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData){
			recalObject.addSite(sites[pos]);
		}
	}
}



