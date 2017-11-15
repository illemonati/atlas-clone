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
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
};

TWindow::TWindow(long Start, long End){
	start = Start;
	end = End;
	initSites(end - start); //end NOT in window!
};


void TWindow::clear(){
	for(int i=0; i<length; ++i) sites[i].clear();
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
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

bool TWindow::addFromRead(BamTools::BamAlignment & bamAlignment, TPMD* pmdObjects, TReadGroups* readGroups, int & minQuality, int & maxQuality){
	/* Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window
	 */
	if(bamAlignment.Position >= end) return true;

	//find first position to be within window
	double len = bamAlignment.AlignedBases.length();
	if(bamAlignment.Position + len < start) return false;

	//find which position to consider first
	++numReadsInWindow;
	int firstPos = start - bamAlignment.Position;
	if(firstPos < 0) firstPos = 0;
	int lastPos = len;
	if(bamAlignment.Position + lastPos > end) lastPos = end - bamAlignment.Position;

	//find relevant 3' end
	int distFrom5Prime, distFrom3Prime;

	//Extract Read Group Info
	std::string readGroup;
	bamAlignment.GetTag("RG", readGroup);
	int readGroupId = readGroups->find(readGroup);

	//add sites
	int internalPos = bamAlignment.Position + firstPos - start;
	char base; BaseContext context;
	char quality;
	int secondLastPos = lastPos - 1;
	double pmdCT, pmdGA;

	/* Note:
	 *  1) Reference is 5' -> 3'
	 *  2) distance is 1-based!
	 *  3) Ignoring indels when calculating distances
	 *  4) Function add needs first P(C->T), then P(G->A)
	 */

	if(bamAlignment.IsProperPair()){
		if(bamAlignment.IsReverseStrand()){
			//reverse (can be either first or second mate, but it's the one that comes second in bam file)
			//hence P(C->T) is given by f(dist since beginning of fragment) = f(insert - len + pos)
			//and P(G->A) is given as f(end of fragment) = f(len - pos - 1)
			for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
				base = bamAlignment.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
					quality = bamAlignment.AlignedQualities.at(pos);
					if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
						//get context: flip bases!
						if(pos == secondLastPos) context = genoMap.getContextReverseRead('N', base);
						else context = genoMap.getContextReverseRead(bamAlignment.AlignedBases.at(pos + 1), base);
						//set distances
						distFrom3Prime = abs(bamAlignment.InsertSize) - len + pos;
						distFrom5Prime = len - pos - 1;
						pmdCT = pmdObjects[readGroupId].getProbGA(distFrom3Prime); //get flipped pmd pattern instead of flipping bases
						pmdGA = pmdObjects[readGroupId].getProbCT(distFrom5Prime);
						//add base
						sites[internalPos].add(base, quality, distFrom5Prime, distFrom3Prime, pmdCT, pmdGA, context, readGroupId);
					}
				}
			}
		} else {
			//forward (can be either first or second mate, but it's the one that comes first in bam file)
			//Hence P(C->T) is given as a function of pos
			//And P(G->A) is given by (length of fragment) - pos -1
			for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
				base = bamAlignment.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
					quality = bamAlignment.AlignedQualities.at(pos);
					if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
						//get context
						if(pos == 0) context = genoMap.getContext('N', base);
						else context = genoMap.getContext(bamAlignment.AlignedBases.at(pos - 1), base);
						//set distances
						distFrom5Prime = pos;
						distFrom3Prime = bamAlignment.InsertSize - pos - 1;
						//add base
						sites[internalPos].add(base, quality, distFrom5Prime, distFrom3Prime, pmdObjects[readGroupId].getProbCT(distFrom5Prime), pmdObjects[readGroupId].getProbGA(distFrom3Prime), context, readGroupId);
					}
				}
			}
		}
	} else {
		//treat as single end
		if(bamAlignment.IsReverseStrand()){
			//not in pair & reverse
			//Hence P(C->T) from 5' is just as P(G->A) from 3' in forward: f(pos)
			//And P(G->A) from 3' is just as P(C->T) from 5' in forward: f(len - pos - 1)
			for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
				base = bamAlignment.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
					quality = bamAlignment.AlignedQualities.at(pos);
					if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
						//get context: flip bases!
						if(pos == secondLastPos) context = genoMap.getContextReverseRead('N', base);
						else context = genoMap.getContextReverseRead(bamAlignment.AlignedBases.at(pos + 1), base);

						//set distances
						distFrom5Prime = len - pos - 1;
						distFrom3Prime = pos;
						//add base
						sites[internalPos].add(base, quality, distFrom5Prime, distFrom3Prime, pmdObjects[readGroupId].getProbGA(distFrom3Prime), pmdObjects[readGroupId].getProbCT(distFrom5Prime), context, readGroupId);
					}
				}
			}
		} else {
			//not in pair & forward
			//Hence P(C->T) is given as a function of pos
			//And P(G->A) is given by len - pos -1
			for(int pos = firstPos; pos < lastPos; ++pos, ++internalPos){
				base = bamAlignment.AlignedBases.at(pos);
				if(base == 'A' || base == 'C' || base == 'G' || base == 'T'){ //skip any other
					quality = bamAlignment.AlignedQualities.at(pos);
					if(minQuality <= (int) quality && (int) quality <= maxQuality){ //skip if quality does not make sense
						//get context
						if(pos == 0) context = genoMap.getContext('N', base);
						else context = genoMap.getContext(bamAlignment.AlignedBases.at(pos - 1), base);
						//set distances
						distFrom5Prime = pos;
						distFrom3Prime = len - pos - 1;
						//add base
						sites[internalPos].add(base, quality, distFrom5Prime, distFrom3Prime, pmdObjects[readGroupId].getProbCT(distFrom5Prime), pmdObjects[readGroupId].getProbGA(distFrom3Prime), context, readGroupId);
					}
				}
			}
		}
	}

	//return if part of the read maps to next window
	if(lastPos == len) return false;
	else return true;
}


bool TWindow::addFromRead(TAlignmentParser & alignemntParser, TPMD* pmdObjects, TReadGroups* readGroups, int & minQuality, int & maxQuality){
	/* Note:
	 * Function returns true if read also maps to next window and
	 * returns false if end of read is within this (or a previous) window
	 */

	//check if alignment is inside window
	if(alignemntParser.position >= end) return true;
	if(alignemntParser.position + alignemntParser.length < start) return false;

	//find which position to consider first
	++numReadsInWindow;
	int firstPos = alignemntParser.position - start;
	int p = 0;
	if(firstPos < 0){
		firstPos = 0;
		while(p < alignemntParser.length && (firstPos + alignemntParser.alignedPos[p]) > 0)
			++p;
		if(p == alignemntParser.length)
			return false;
	}
	int internalPos;

	/* Note:
	 *  1) Reference is 5' -> 3'
	 *  2) distance is 0-based!
	 *  3) Ignoring indels in other mate when calculating distances
	 *  4) Function add needs first P(C->T), then P(G->A)
	 */

	for(; p < alignemntParser.length; ++p){
		if(alignemntParser.aligned[p] && alignemntParser.base[p] != N){
			if(minQuality <= alignemntParser.quality[p] && alignemntParser.quality[p] <= maxQuality){ //skip if quality does not make sense
				internalPos = firstPos + alignemntParser.alignedPos[p];
				if(internalPos >= length)
					return true; //since part of the read maps to next window
				sites[internalPos].add(alignemntParser.base[p], alignemntParser.quality[p], alignemntParser.distFrom5Prime[p], alignemntParser.distFrom3Prime[p], pmdObjects[alignemntParser.readGroupId].getProbCT(alignemntParser.distFrom5Prime[p]), pmdObjects[alignemntParser.readGroupId].getProbGA(alignemntParser.distFrom3Prime[p]), alignemntParser.context[p], alignemntParser.readGroupId);
			}
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
		if(sites[i].hasData){
			sites[i].addToBaseFrequencies(baseFreq);
		}
	}
	baseFreq.normalize();
}

void TWindow::calculateEmissionProbabilities(TRecalibration* recalObject){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject->calcEmissionProbabilities(sites[i]);
		}
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

void TWindow::printPileup(TRecalibration* recalObject, std::ofstream & out, std::string & chr){
	//calc emission probs
	for(int i=0; i<length; ++i){
		recalObject->calcEmissionProbabilities(sites[i]);
	}
	//print pileup
	for(int i=0; i<length; ++i){
		out << chr << "\t" << start + i + 1 << "\t" << sites[i].bases.size() << "\t" << sites[i].getBases() << "\t" << sites[i].getEmissionProbs() << "\n";
	}
}

void TWindow::calcCoverage(){
	//calculate and return coverage
	coverage = 0.0;
	long noData = 0;
	long plentyData = 0;
	int cov;
	for(int i=0; i<length; ++i){
		cov = sites[i].depth();
		coverage += cov;
		if(cov == 0) ++noData;
		else if(cov > 1) ++ plentyData;
	}

	coverage = coverage / (double) length;
	numSitesWithData = length - noData;
	fractionSitesNoData = (double) noData / (double) length;
	fractionsitesCoverageAtLeastTwo = (double) plentyData / (double) length;
}

void TWindow::calcFracN(){
	double numN = 0.0;
	for(int i=0; i<length; ++i)	if(sites[i].referenceBase == 'N') ++numN;
	fractionRefIsN = numN / (double) length;
}

void TWindow::calcCoveragePerSite(long* siteCoverage, size_t maxCov){
	//calculate and return coverage
	coverage = 0.0;
	long noData = 0;
	long plentyData = 0;
	size_t cov;

	for(int i=0; i<length; ++i){
		cov = sites[i].depth();
		if(cov <= maxCov)	siteCoverage[cov] += 1;
		else siteCoverage[maxCov + 1] += 1; //else it should be in the "greater than" bin

		if(cov == 0) ++ noData;
		else if(cov > 1) ++ plentyData;
	}

	coverage = coverage / (double) length;
	fractionSitesNoData = (double) noData / (double) length;
	fractionsitesCoverageAtLeastTwo = (double) plentyData / (double) length;
}

void TWindow::applyCoverageFilter(int minCoverage, size_t maxCoverage){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minCoverage || sites[i].bases.size() > maxCoverage)
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
				QTtables.at((*it)->readGroup)->add((*it)->quality, recalObject->getQuality(*it));
				QTtables.at(QTtables.size() - 1)->add((*it)->quality, recalObject->getQuality(*it));
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
				QTtables.at((*it)->readGroup)->add(recalObject->getQuality(*it), otherRecalObject->getQuality(*it));
				QTtables.at(QTtables.size() - 1)->add(recalObject->getQuality(*it), otherRecalObject->getQuality(*it));
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
void TWindowDiploid::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	try{
		sites = new TSiteDiploid[length];
	} catch(...){
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
	}
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
}

void TWindowDiploid::addSitesToThetaEstimator(TRecalibration* recalObject, TThetaEstimator & estimator){
	//first calculate emission probabilities
	calculateEmissionProbabilities(recalObject);

	//now add sites
	addSitesToThetaEstimator(estimator);
}

void TWindowDiploid::addSitesToThetaEstimator(TThetaEstimator & estimator){
	//assumes that emission probabilities were calculated
	for(int i=0; i<length; ++i){
		estimator.add(sites[i]);
	}
}

void TWindowDiploid::callMLEGenotypeKnownAlleles(TRecalibration* recalObject, TSiteSubset* subset, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool & isVCF, bool & noAltIfHomoRef, bool & beagle, bool & printOnlyGL){
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

void TWindowDiploid::callBayesianGenotype(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF){
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

void TWindowDiploid::callBayesianGenotypeKnownAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF){
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

void TWindowDiploid::callAllelePresence(TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool noAltIfHomoRef){
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

void TWindowDiploid::callRandomBase(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll){
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

void TWindowDiploid::majorityCall(TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll){
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

void TWindowDiploid::callAllelePresenceKnwonAlleles(TSiteSubset* subset, TThetaEstimator & estimator, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool isVCF, bool noAltIfHomoRef){
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

void TWindowDiploid::addToGLF(TGlfWriter & writer, bool printAll){
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

void TWindowDiploid::generatePSMCInput(TThetaEstimator & estimator, int & blockSize, double & confidence, std::ofstream & out, int & nCharOnLine){
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

//-------------------------------------------------------
//TWindowHaploid
//-------------------------------------------------------

void TWindowHaploid::initSites(long newLength){
	if(sitesInitialized)
		delete[] sites;
	length = newLength;
	try{
		sites = new TSiteHaploid[length];
	} catch(...){
		throw "Failed to allocate sufficient memory to store the data for so many sites. Consider reducing the window size or selecting fewer sites.";
	}
	sitesInitialized = true;
	coverage = -1.0;
	fractionSitesNoData = -1.0;
	fractionsitesCoverageAtLeastTwo = -1.0;
	numReadsInWindow = 0;
}

void TWindowHaploid::fillPGenotype(double* pGenotype){
	for(int i=0; i<4; ++i){
		pGenotype[i] = baseFreq[i];
	}
}

double TWindowHaploid::calcLogLikelihood(){
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

void TWindowHaploid::addToRecalibrationEM(TRecalibrationEM & recalObject){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject.addSite(sites[i]);
		}
	}
}

void TWindowHaploid::addToRecalibrationEM(TRecalibrationEM & recalObject, TSiteSubset* subset){
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

void TWindowHaploid::addToExpectedBaseCounts(TRecalibration* recalObject, double** expectedCounts){
	estimateBaseFrequencies();
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject->calcEmissionProbabilities(sites[i]);
			sites[i].addToExpectedBaseCounts(baseFreq, expectedCounts[i]);
		}
	}
}

void TWindowHaploid::calculatePoolFreqLikelihoods(int & numChromosomes, Base** majorMinor, gz::ogzstream & out, std::string & chr, bool printAll){
	//assumes that emission probabilities were calculated!!
	if(printAll){
		for(int i=0; i<length; ++i){
			out << chr << "\t" << start + i + 1;
			sites[i].calculatePoolFreqLikelihoods(numChromosomes, genoMap, majorMinor[i][0], majorMinor[i][1], out);
		}
	} else {
		for(int i=0; i<length; ++i){
			if(sites[i].hasData){
				out << chr << "\t" << start + i + 1;
				sites[i].calculatePoolFreqLikelihoods(numChromosomes, genoMap, majorMinor[i][0], majorMinor[i][1], out);
			}
		}
	}
}

