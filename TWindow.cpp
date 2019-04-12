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
//	chrIterator = NULL;
	chrNumber = -1;
	sites = NULL;
	sitesInitialized = false;
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesDepthAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
	//lastAlignmentwithEndInWindow = usedAlignments.end();
	//firstAlignmentwithPosOutsideWindow = usedAlignments.end();
	passedFilters = false;
};

TWindow::~TWindow(){
	//delete sites
	clear();
	if(sitesInitialized)
		delete[] sites;

	//delete alignments
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	usedAlignments.clear();

	for(std::vector<TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		delete *alignmentIt;
	emptyAlignments.clear();

};

TAlignment* TWindow::swapUsedForEmptyAlignment(TAlignment* usedAlignment, const unsigned int & maxReadLength){
	//save used alignment on proper stack
	usedAlignments.push_back(usedAlignment);

	//return empty alignment, either from stack or create new
	if(emptyAlignments.size() > 0){
		TAlignment* alignment = *(emptyAlignments.rbegin());
		emptyAlignments.pop_back();
		return alignment;
	} else {
		TAlignment* alignment = new TAlignment(maxReadLength);
		return alignment;
	}
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
		for(int i=0; i<length; ++i)
			sites[i].clear();
	}
	depth = -1.0;
	fractionSitesNoData = -1.0;
	fractionRefIsN = -1.0;
	fractionsitesDepthAtLeastTwo = -1.0;
	numSitesWithData = 0;
	numReadsInWindow = 0;
	referenceBaseAdded = false;
	passedFilters = false;
};

void TWindow::setCoordinates(long Start, long End, int ChrNumber){
	start = Start;
	end = End;
	chrNumber = ChrNumber;
	if(sitesInitialized){
		if((end - start) != length)
			initSites(end - start);
		else
			clear();
	} else initSites(end - start);
}

void TWindow::move(long Start, long End, int ChrNumber){
	setCoordinates(Start, End, ChrNumber);
	cleanUpUsedAlignments();
};

void TWindow::jump(long Start, long End, int ChrNumber){
	setCoordinates(Start, End, ChrNumber);
	clearAllUsedAlignments();
}

void TWindow::review(){
	//update pointers
	/*
	firstAlignmentwithPosOutsideWindow = usedAlignments.end()-1;
	while((*firstAlignmentwithPosOutsideWindow)->position > end && firstAlignmentwithPosOutsideWindow != usedAlignments.begin())
		--firstAlignmentwithPosOutsideWindow;
		*/

	//fillSites();
	//calcDepth();
}

void TWindow::cleanUpUsedAlignments(){
//	std::cout << "cleaning up used alignemtns:" << std::flush;
	//now check and move the rest
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		if((*alignmentIt)->position < end && (*alignmentIt)->lastAlignedPositionWithRespectToRef >= start && (*alignmentIt)->chrNumber == chrNumber){
			++alignmentIt;
		} else{
			(*alignmentIt)->clear();
			emptyAlignments.push_back(*alignmentIt);
			alignmentIt = usedAlignments.erase(alignmentIt);
		}

//		if((*alignmentIt)->position >= end){
//			usedAlignments.erase(alignmentIt, usedAlignments.end());
////			std::cout << (*alignmentIt)->alignmentName << ":" << (*alignmentIt)->position<< "\t" << std::flush;
//			break;
//		}
//		std::cout << std::endl;
	}
}

void TWindow::clearAllUsedAlignments(){
	for(std::vector<TAlignment*>::iterator alignmentIt=usedAlignments.begin(); alignmentIt != usedAlignments.end();){
		(*alignmentIt)->clear();
		emptyAlignments.push_back(*alignmentIt);
		alignmentIt = usedAlignments.erase(alignmentIt);
	}
}

void TWindow::printStacks(){
	std::cout << "USED ALIGMENTS:";
	for(TAlignment* alignmentIt : usedAlignments)
		std::cout << " " << alignmentIt << " : " << alignmentIt->alignmentName << " pos " << alignmentIt->position;
	std::cout << std::endl;

	std::cout << "EMPTY ALIGMENTS:";
	for(std::vector<TAlignment*>::iterator alignmentIt=emptyAlignments.begin(); alignmentIt != emptyAlignments.end(); ++alignmentIt)
		std::cout << " " << *alignmentIt;
	std::cout << std::endl;

}

void TWindow::fillSitesSubset(TSiteSubset* subset, const int & readUpToDepth){
	//add reads in usedAlignments to sites in window
	for(TAlignment* alignmentIt : usedAlignments){
		//check if alignment start is inside window
		if(alignmentIt->position >= end)
			throw "alignment should be assigned to next window!";

		//genomic position of alignment as seen from window perspective
		int firstPos = alignmentIt->position - start;

		//set position in read
		int p = 0;

		//is the beginning of the read part of previous window? increase starting p for adding bases!
		if(firstPos < 0){
			while(p < alignmentIt->length && (firstPos + alignmentIt->bases[p].alignedPos) < 0)
				++p;
			if(p == alignmentIt->length){
				throw "alignment should be assigned to previous window! Name: " + alignmentIt->alignmentName + ". In window " + toString(start) + "-" + toString(end) + ". with position " + toString(alignmentIt->position);
			}
		}

		//get positions that are used
		std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);

		//position in window where first one = 0
		int internalPos;
		//p is at first position of read in window
		for(; p < alignmentIt->length; ++p){
			if(alignmentIt->bases[p].alignedPos && alignmentIt->bases[p].base != N){
				internalPos = firstPos + alignmentIt->bases[p].alignedPos;
				//if read extends past window length
				if(internalPos >= length)
					break; //since part of the read maps to next window
				if(thesePos.find(internalPos) != thesePos.end() && sites[internalPos].depth() < readUpToDepth)
					sites[internalPos].add(&alignmentIt->bases[p]);
			}
		}
		++numReadsInWindow;
	}
};

void TWindow::fillSites(const int & readUpToDepth){
	//add reads in usedAlignments to sites in window
	for(TAlignment* alignmentIt : usedAlignments){
		//check if alignment start is inside window
		if(alignmentIt->position >= end){
			throw "alignment should be assigned to next window!";
		}

		//genomic position of alignment as seen from window perspective
		int firstPos = alignmentIt->position - start;

		//set position in read
		int p = 0;

		//is the beginning of the read part of previous window? increase starting p for adding bases!
		if(firstPos < 0){
			while(firstPos + alignmentIt->bases[p].alignedPos < 0){
				++p;
				if(p == alignmentIt->length){
					throw "alignment should be assigned to previous window! Name: " + alignmentIt->alignmentName + ". In window " + toString(start) + "-" + toString(end) + ". with position " + toString(alignmentIt->position);
				}
			}
		}

		//position in window where first one = 0
		int internalPos;
		//p is at first position of read in window
		for(; p < alignmentIt->length; ++p){
				if(alignmentIt->bases[p].aligned && alignmentIt->bases[p].base != N){
					internalPos = firstPos + alignmentIt->bases[p].alignedPos;
					//if read extends past window length
					if(internalPos >= length)
						break; //since part of the read maps to next window
					if(sites[internalPos].depth() < readUpToDepth)
						sites[internalPos].add(&alignmentIt->bases[p]);
				}
		}
		++numReadsInWindow;
	}
};

void TWindow::addReferenceBaseToSites(BamTools::Fasta & reference, int & refId){
	if(!referenceBaseAdded){
		int stop = end - 1; //note that end is last position + 1
		std::string ref; //fasta object fills string
		reference.GetSequence(refId, start, stop, ref);
		for(int i=0; i<length; ++i){
			sites[i].setRefBase(ref[i]);
		}
		referenceBaseAdded = true;
	}
}

void TWindow::addReferenceBaseToSites(TSiteSubset* subset){
	if(!referenceBaseAdded){
		std::cout << "ADDING REF TO SITES!" << std::endl;
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

		std::cout << "SUCCESSFULL!" << std::endl;
	}
}

void TWindow::applyMask(TBedReader* mask, bool doInverseMasking){
	int pos;
	long first = 0;
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

void TWindow::calculateEmissionProbabilities(){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData)
			sites[i].calcEmissionProbabilities();
	}
}

void TWindow::call(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference){
	//add reference to sites
	addReferenceBaseToSites(reference, chrNumber);

	//loop over sites and call
	for(int i=0; i<length; ++i){
		if(sites[i].hasData)
			sites[i].calcEmissionProbabilities();
		caller.call(chrName, start + i + 1, sites[i]); //i + 1 to make vcf 1-based!
	}
};

void TWindow::callKnwonAlleles(TCaller & caller, TRecalibration & recalObject, BamTools::Fasta & reference, TSiteSubset & subset){
	//check if we need to process this window
	if(subset.hasPositionsInWindow(start)){
		//add reference to sites
		addReferenceBaseToSites(&subset);

		//now only run over sites listed in that window
		std::map<long,std::pair<char,char> > thesePos = subset.getPositionInWindow(start);
		for(std::map<long,std::pair<char,char> >::iterator it = thesePos.begin(); it!=thesePos.end(); ++it){
			int pos = it->first - start;
			if(sites[pos].hasData){
				sites[pos].calcEmissionProbabilities();

				//call
				caller.call(chrName, start + pos + 1, sites[pos], it->second.first, it->second.second); //pos + 1 to make vcf 1-based
			}
		}
	}
};



/*

void TWindow::callMLEGenotype(TRecalibration* recalObject, TRandomGenerator & randomGenerator, gz::ogzstream & out, std::string & chr, bool printAll, bool printRef, bool isVCF, bool gVCF, bool noAltIfHomoRef){
	if(isVCF){
		if(printAll){
			for(int i=0; i<length; ++i){
				out << chr << "\t" << start + i + 1;
				if(sites[i].hasData)
					sites[i].calcEmissionProbabilities();
				std::string basesString = sites[i].getBases();
				sites[i].callMLEGenotypeVCF(genoMap, randomGenerator, out, gVCF, noAltIfHomoRef, basesString);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					sites[i].calcEmissionProbabilities();
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
				if(sites[i].hasData)
					sites[i].calcEmissionProbabilities();
				sites[i].callMLEGenotype(genoMap, randomGenerator, out);
				out << "\n";
			}
		} else {
			for(int i=0; i<length; ++i){
				if(sites[i].hasData){
					out << chr << "\t" << start + i + 1;
					sites[i].calcEmissionProbabilities();
					sites[i].callMLEGenotype(genoMap, randomGenerator, out);
					out << "\n";
				}
			}
		}
	}
}

*/
void TWindow::printPileup(TRecalibration* recalObject, gz::ogzstream & out, bool printOnlySitesWithData){
	//print pileup
	for(int i=0; i<length; ++i){
		if((printOnlySitesWithData && sites[i].hasData) || !printOnlySitesWithData){
			sites[i].calcEmissionProbabilities();
			out << chrName << "\t" << start + i + 1;
			sites[i].printPileup(out);
			out << "\n";
		}
	}
}

void TWindow::printPileupToScreen(TRecalibration* recalObject){
	//print pileup
	for(int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		std::cout << chrName << "\t" << start + i + 1;
		sites[i].printPileupToScreen();
		std::cout << "\n";
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
};

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
};

void TWindow::countDepthPerSite(TDistributionOfCounts & counts){
	for(int i=0; i<length; ++i)
		counts.add(sites[i].depth());
};

void TWindow::printDepthPerSite(gz::ogzstream & out, std::string & chr){
	//print depth for each site to file
	for(int i=0; i<length; ++i){
		out << chr << "\t" << start + i + 1 << "\t" << sites[i].depth() << "\n";
	}
};

void TWindow::printMateInformationPerSite(TOutputFileZipped & out, std::string & chr){
	int* alleleCounts = new int[4];
	int* mateCounts = new int[2];
	int* frCounts = new int[2];

	for(int i=0; i<length; ++i){
		//chr, pos and depth
		out << chr << start + i + 1 << sites[i].depth();

		//num alleles
		sites[i].countAlleles(alleleCounts);
		int numAlleles = 0;
		for(int b=0; b<4; ++b){
			out << alleleCounts[b];
			if(alleleCounts[b] > 0)
				++numAlleles;
		}
		out << numAlleles;

		//mates
		sites[i].countMates(mateCounts);
		out << mateCounts[0] << mateCounts[1];

		//fwd / rev
		sites[i].countFwdRev(frCounts);
		out << frCounts[0] << frCounts[1];
		out << std::endl;
	}

	delete[] alleleCounts;
	delete[] mateCounts;
	delete[] frCounts;
};

void TWindow::countAlleles(long**** siteImbalance, const unsigned int & maxCov){
	//calculate and return imbalance
	for(int i=0; i<length; ++i){
		if(sites[i].depth() <= maxCov && sites[i].depth() > 0)
			sites[i].countAllelesForImbalance(siteImbalance);
		else if(sites[i].depth()  == 0){
			++siteImbalance[0][0][0][0];
		}

	}
};

void TWindow::applyDepthFilter(const size_t minDepth, const size_t maxDepth){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepth || sites[i].bases.size() > maxDepth)
				sites[i].clear();
		}
	}
};

void TWindow::createDepthMask(size_t minDepthForMask, size_t maxDepthForMask, std::ofstream & outputMaskFile, std::string & chr){
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			if(sites[i].bases.size() < minDepthForMask || sites[i].bases.size() > maxDepthForMask){
				outputMaskFile << chr << "\t" << start + i << "\t" << start + i + 1 << "\n";
			}
		}
	}
};

void TWindow::addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TLog* logfile, TQualityMap & qualMap){
	logfile->listFlush("Adding sites to BQSR ...");
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			bqsr.addSite(sites[i], qualMap);
		}
	}
	logfile->done();
};

void TWindow::addSitesToBQSR(TRecalibrationBQSREstimator & bqsr, TSiteSubset* subset, TLog* logfile, TQualityMap & qualMap){
	logfile->listFlush("Adding sites to BQSR ...");
	//now only run over sites listed in that window
	std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData){
			sites[pos].setRefBase(it->second.second);
			bqsr.addSite(sites[pos], qualMap);
		}
	}
	logfile->done();

};


void TWindow::addSitesToPMDTable(TPMDTables & pmdTables, TLog* logfile){
	logfile->listFlush("Adding sites to PMD tables ...");
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			//pmdTables.add(sites[i]);
		}
	}
	logfile->done();
};

void TWindow::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer){
	//assumes that emission probabilities were calculated
	for(int i=0; i<length; ++i){
		sites[i].calcEmissionProbabilities();
		thetaDataContainer->add(sites[i]);
	}
};

void TWindow::addSitesToThetaEstimator(TThetaEstimatorData* thetaDataContainer, TBedReader & region){
	//assumes that emission probabilities were calculated
	//only add sites from regions
	if(region.hasPositionsInWindow(start)){
		std::vector<long> thesePos = region.getPositionInWindow(start);
		for(std::vector<long>::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
			int pos = *it - start;
			if(pos < length){
				sites[pos].calcEmissionProbabilities();
				thetaDataContainer->add(sites[pos]);
			}
		}
	}
};

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
};

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

void TWindow::addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	for(int i=0; i<length; ++i){
		if(sites[i].hasData){
			recalObject.addSite(sites[i], qualMap);
		}
	}
}

void TWindow::addToRecalibrationEM(TRecalibrationEMEstimator & recalObject, TSiteSubset* subset, TQualityMap & qualMap){
	estimateBaseFrequencies();
	recalObject.addNewWindow(&baseFreq);
	//now only run over sites listed in that window
	std::map<long,std::pair<char,char> > thesePos = subset->getPositionInWindow(start);
	int pos;
	for(std::map<long,std::pair<char,char> >::iterator it=thesePos.begin(); it!=thesePos.end(); ++it){
		pos = it->first - start;
		if(sites[pos].hasData){
			recalObject.addSite(sites[pos], qualMap);
		}
	}
}



