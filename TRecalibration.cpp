/*
 * TRecalibration.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#include "TRecalibration.h"

//---------------------------------------------------------------
//RecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(TParameters* arguments, TLog* logfile){
	numParams = 3; //3 beta and 4 gamma
	params = new double[numParams];
	newParams = new double[numParams];
	Jacobian.resize(numParams,numParams);
	Jacobian.zeros();
	F.resize(numParams);
	F.zeros();
	JxF.resize(numParams, numParams);
	JxF.zeros();
	maxF = 0.0;
	numSitesAdded = 0;
	logLikelihood = 0.0;

	//set initial values
	//betas
	logfile->startIndent("Will start recalibration EM with:");
	params[0] = arguments->getParameterDoubleWithDefault("initBeta0", 0);
	logfile->list("beta0 = " + toString(params[0]));
	params[1] = arguments->getParameterDoubleWithDefault("initBeta1", 1);
	logfile->list("beta1 = " + toString(params[1]));
	params[2] = arguments->getParameterDoubleWithDefault("initBeta2", 0);
	logfile->list("beta2 = " + toString(params[2]));

	//gammas -> start all at 0 by default
	/*
	params[3] = arguments->getParameterDoubleWithDefault("initGammaA", 0);
	logfile->list("gammaA = " + toString(params[3]));
	params[4] = arguments->getParameterDoubleWithDefault("initGammaC", 0);
	logfile->list("gammaC = " + toString(params[4]));
	params[5] = arguments->getParameterDoubleWithDefault("initGammaG", 0);
	logfile->list("gammaG = " + toString(params[5]));
	params[6] = arguments->getParameterDoubleWithDefault("initGammaT", 0);
	logfile->list("gammaT = " + toString(params[6]));
	*/
	logfile->endIndent();
};

TRecalibrationEM::TRecalibrationEM(TLog* logfile){
	numParams = 3; //3 beta and 4 gamma
	params = new double[numParams];
	newParams = new double[numParams];
	Jacobian.resize(numParams,numParams);
	Jacobian.zeros();
	F.resize(numParams);
	F.zeros();
	JxF.resize(numParams, numParams);
	JxF.zeros();
	maxF = 0.0;
	numSitesAdded = 0;
	logLikelihood = 0.0;

	//init params
	for(int i=0; i<numParams; ++i){
		params[i] = 0.0;
		newParams[i] = 0.0;
	}
}

void TRecalibrationEM::setParams(double* Params){
	for(int i=0; i<numParams; ++i){
		params[i] = Params[i];
	}
}

double TRecalibrationEM::calcEta(TBase* base){
	return calcEta(base, params);
}

double TRecalibrationEM::calcEta(TBase* base, double* theseParams){
	//function transform log error
	double eta = theseParams[0] + theseParams[1] * base->transformedLogError + theseParams[2] * (double) base->posInRead;
	//eta += theseParams[base->getBaseAsEnum() + 3];
	return eta;
}

double TRecalibrationEM::calcEpsilon(const double & eta){
	double tmp = exp(eta);
	return tmp / (1.0 + tmp);
}

double TRecalibrationEM::calcEpsilon(TBase* base){
	return calcEpsilon(base, params);
}

double TRecalibrationEM::calcEpsilon(TBase* base, double* theseParams){
	double tmp = calcEta(base, theseParams);
	return calcEpsilon(tmp);
}

void TRecalibrationEM::initEMStep(){
	//fill vector of new params by copying current values
	for(int i=0; i<numParams; ++i){
		newParams[i] = params[i];
	}
}

void TRecalibrationEM::initNetwonRalphsonStep(){
	Jacobian.zeros();
	F.zeros();
	numSitesAdded = 0;
}

void TRecalibrationEM::saveParams(){
	for(int i=0; i<numParams; ++i){
		params[i] = newParams[i];
	}
}

void TRecalibrationEM::addSiteToJacobianAndF(std::vector<TBase*> & bases, TBaseFrequencies* freqs){
	//adds terms to Jacobian for one site (hence a vector of bases that were read)
	//assumes bases to be haploid! -> program does not throw an error if they are not!
	double epsilon;
	double epsilonThird;

	//initialize
	double P_d_given_g_theta[4];
	for(int g=0; g<4; ++g){ //over all genotypes A, C, G or T
		P_d_given_g_theta[g] = 1.0;
	}

	//calc P(d|g, theta) for all g
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		//get epsilon
		epsilon  = calcEpsilon(*it); //use current (= old) params for this step
		epsilonThird = epsilon / 3.0;

		for(int g=0; g<4; ++g){ //over all genotypes
			//add to P(d|g, theta)
			if((*it)->getBaseAsEnum() == g)
				P_d_given_g_theta[g] *= (1.0 - epsilon);
			else
				P_d_given_g_theta[g] *= epsilonThird;
		}
	}

	//calculate P(g|d, theta)
	double P_g_given_d_theta_denominator = 0.0;
	for(int g=0; g<4; ++g){
		P_d_given_g_theta[g] *= (*freqs)[g];
		P_g_given_d_theta_denominator += P_d_given_g_theta[g];
	}
	double P_g_given_d_theta[4];
	for(int g=0; g<4; ++g){
		P_g_given_d_theta[g] = P_d_given_g_theta[g] / P_g_given_d_theta_denominator;
	}

	//now add site to Jacobian and F
	//Note: Jacobian is symmetric! But we only add to top triangle -> will copy final numbers down later
	int g, nucIndex;
	double eta, expEta, tmp, onePlusExpEta, weightFcorrect, weightFincorrect;
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		eta = calcEta(*it, newParams); //use new params for this step!
		expEta = exp(eta);
		onePlusExpEta = 1.0 + expEta;
		tmp = expEta / (onePlusExpEta * onePlusExpEta);

		Jacobian(0,0) -= tmp;
		Jacobian(0,1) -= tmp * (*it)->transformedLogError;
		Jacobian(0,2) -= tmp * (double) (*it)->posInRead;
		Jacobian(1,1) -= tmp * (*it)->transformedLogError * (*it)->transformedLogError;
		Jacobian(1,2) -= tmp * (*it)->transformedLogError * (double) (*it)->posInRead;
		Jacobian(2,2) -= tmp * (double) (*it)->posInRead * (double) (*it)->posInRead;

		//now add to those with indicators on the base
		//Note that the derivative is zero for all other bases
		g = (*it)->getBaseAsEnum();
		nucIndex = g + 3;
		/*
		Jacobian(nucIndex, nucIndex) -= tmp;
		Jacobian(0, nucIndex) -= tmp;
		Jacobian(1, nucIndex) -= tmp * (*it)->transformedLogError;
		Jacobian(2, nucIndex) -= tmp * (*it)->posInRead;
		 */

		//now add to F -> need to add for all genotypes, but with different term
		weightFcorrect = -1.0 / (1.0 + exp(-eta));
		weightFincorrect = 1.0 / onePlusExpEta;
		for(int geno = 0; geno < 4; ++geno){
			//prepare weight
			if(geno == g) tmp = weightFcorrect;
			else tmp = weightFincorrect;
			tmp = P_g_given_d_theta[geno] * tmp;

			//fill for beta0, beta1 and beta2
			F(0) += tmp; //beta0
			F(1) += tmp * (*it)->transformedLogError;
			F(2) += tmp * (double) (*it)->posInRead;
		}
		//derivation is only non-zero for beta3 to beta6 if g corresponds to index
		//F(nucIndex) += P_g_given_d_theta[g] * weightFcorrect;
	}

	++numSitesAdded;
}

void TRecalibrationEM::runNewtonRalphson(){
	//Need to copy numbers to other triangle in Jacobian, as only upper triangled is filled when parsing sites
	for(int i=0; i<(numParams-1); ++i){
		for(int j=i+1; j<numParams; ++j){
			//copy from upper triangle to lower triangle
			Jacobian(j,i) = Jacobian(i,j);
		}
	}

	//scale F and J by 1/#sites
	Jacobian = Jacobian / (double) numSitesAdded;
	F = F / (double) numSitesAdded;

	std::cout << "JACOBIAN = " << std::endl;
	std::cout << Jacobian << std::endl;

	std::cout << std::endl << "DET = " << det(Jacobian) << std::endl;

	std::cout << "F = " << std::endl;
		std::cout << F << std::endl;

	//get largest gradient (F) to check if we break
	maxF = 0.0;
	for(int i=0; i<(numParams-1); ++i){
		if(F(i) > maxF) maxF = F(i);
	}
	std::cout << "MAX F = " << maxF << std::endl;

	//now calculate J^-1 x F
	std::cout << "New Params:";

	if(solve(JxF, Jacobian, F)){

		std::cout << "JxF = " << std::endl;
		std::cout << JxF << std::endl;

		for(int i=0; i<numParams; ++i){
			newParams[i] -= JxF(i);
			std::cout << "\t" << newParams[i];
		}
	} else throw "Issue solving JxF in TRecalibrationEM::runNewtonRalphson()!";
	std::cout << std::endl;
}

void TRecalibrationEM::writeHeader(std::ofstream & out){
	out << "beta0\tbeta1\tbeta2";
}

void TRecalibrationEM::writeParams(std::ofstream & out){
	out << params[0];
	for(int i=1; i<numParams; ++i){
		out << "\t" << params[i];
	}
}

void TRecalibrationEM::resetLikelihood(){
	logLikelihood = 0.0;
}

void TRecalibrationEM::addSiteToLikelihood(std::vector<TBase*> & bases, TBaseFrequencies* freqs){
	//adds this site (hence a vector of bases that were read) to the log likelihood
	//assumes bases to be haploid! -> program does not throw an error if they are not!
	double epsilon;
	double epsilonThird;
	double* likelihood = new double[4]; //for each genotype
	for(int i=0; i<4; ++i) likelihood[i] = 1.0;

	//calc likelihood term
	for(std::vector<TBase*>::iterator it = bases.begin(); it != bases.end(); ++it){
		//get epsilon
		epsilon  = calcEpsilon(*it); //use current params for this step
		epsilonThird = epsilon / 3.0;

		for(int g=0; g<4; ++g){ //over all genotypes
			if((*it)->getBaseAsEnum() == g){
				likelihood[g] *= (1.0 - epsilon);
			} else {
				likelihood[g] *= epsilonThird;
			}
		}
	}

	//build weighted sum across genotypes
	double ll = 0.0;
	for(int g=0; g<4; ++g){
		ll += (*freqs)[g] * likelihood[g];
	}
	logLikelihood += log(ll);
}


//---------------------------------------------------------------
//TRecalibrationBQSR_cell BQSR
//---------------------------------------------------------------
TBQSR_cellQuality::TBQSR_cellQuality(){
	curEstimate = 0.0;
	estimationConverged = false;
	firstDerivative = 0.0;
	secondDerivative = 0.0;
	pmdObject = NULL;
	numObservations = 0;
	numMatches = 0;
	F = 0.0;
}

void TBQSR_cellQuality::init(double initialError, TPMD* PmdObject){
	curEstimate = initialError;
	pmdObject = PmdObject;
}

void TBQSR_cellQuality::empty(){
	if(!estimationConverged){
		numObservations = 0;
		numMatches = 0;
		firstDerivative = 0.0;
		secondDerivative = 0.0;
	}
}

void TBQSR_cellQuality::addToDerivatives(TBase* base, Base & RefBase, double & epsilon){
	double D = 0.0;
	switch(base->getBaseAsEnum()){
		case A: if(RefBase == A){
					D = 1.0;
					break;
				}
				if(RefBase == G) D = pmdObject->getProbGA(base->pos3);
				break;
		case C: if(RefBase == C) D = pmdObject->getProbGA(base->pos5);
				break;
		case G: if(RefBase == G) D = 1.0 - pmdObject->getProbGA(base->pos3);
				break;
		case T: if(RefBase == C){
					D = pmdObject->getProbGA(base->pos5);
					break;
				}
				if(RefBase == T) D = 1.0;
				break;
		case N: throw "Can not add base with unknown reference to BQSR cell!";
	}

	double oneMinus4D = 1.0 - 4.0 * D;
	firstDerivative += oneMinus4D / (-4.0*D*epsilon + 3.0*D + epsilon);
	double tmp = oneMinus4D / (D*(3.0-4.0*epsilon) + epsilon);
	secondDerivative += tmp * tmp;

	++numObservations;
	if(base->getBaseAsEnum() == RefBase) ++numMatches;
}


void TBQSR_cellQuality::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		addToDerivatives(base, RefBase, curEstimate);
	}
}

bool TBQSR_cellQuality::estimate(double & convergenceThreshold){
	if(estimationConverged) return estimationConverged;

	//check if we only made errors or no errors, in which case epsilon is easily estimated
	if(numMatches == numObservations){
		curEstimate = 0.0;
		estimationConverged = true;
		return estimationConverged;
	}
	if(numMatches == 0){
		curEstimate = 1.0;
		estimationConverged = true;
		return estimationConverged;
	}

	//Else, need Newton-Ralphson to estimate epsilon
	std::cout << "ESTIMATE: " << curEstimate << " + " << firstDerivative / secondDerivative << " (" << firstDerivative << " / " <<  secondDerivative << ") = ";

	curEstimate = curEstimate + firstDerivative / secondDerivative;

	std::cout << curEstimate << std::endl;

	//decide on convergence
	F = fabs(firstDerivative);
	if(F < convergenceThreshold) estimationConverged = true;
	return estimationConverged;
}

//---------------------------------------------------------------
TBQSR_cellPosition::TBQSR_cellPosition():TBQSR_cellQuality(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
}

void TBQSR_cellPosition::init(TBQSR_cellQuality** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject){
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;
	pmdObject = PmdObject;
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate * curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPosition(){
	BQSR_quality_position = NULL;
	considerPosition = false;
}

void TBQSR_cellContext::init(TBQSR_cellQuality** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex, TPMD* PmdObject){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, PmdObject);
	BQSR_quality_position = gotBQSR_quality_position;
	considerPosition = true;
}

void TBQSR_cellContext::init(TBQSR_cellQuality** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex, TPMD* PmdObject){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex, PmdObject);
	BQSR_quality_position = NULL;
}

void TBQSR_cellContext::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		int qIndex = qualityIndex->getIndex(base->quality);
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qIndex].curEstimate * curEstimate;
		if(considerPosition) epsilonAlpha *= BQSR_quality_position[qIndex][base->posInRead].curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TPMD* PmdObject, TLog* Logfile){
	logfile = Logfile;
	bamHeader = BamHeader;
	pmdObject = PmdObject;
	numReadGroups = bamHeader->ReadGroups.Size();
	estimationConverged = false;
	qualityConverged = false;
	positionConverged = false;
	contextConverged = false;
	numContexts = 25;

	//read arguments from user
	int minQ = params.getParameterIntWithDefault("minQ", 1);
	int maxQ = params.getParameterIntWithDefault("maxQ", 100);
	logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
	qualityIndex = new TQualityIndex(minQ, maxQ);
	initialError = params.getParameterDoubleWithDefault("initialError", 0.01);
	logfile->list("Using an initial error of " + toString(initialError));
	convergenceThreshold = params.getParameterDoubleWithDefault("maxF", 0.0001);
	logfile->list("Stopping Newton-Ralphson if F < " + toString(convergenceThreshold));

	//read which covariates to consider
	considerPosition = params.parameterExists("considerPosition");

	//initialize BQSR table
	BQSR_cells_readGroup_quality = new TBQSR_cellQuality*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_quality[i] = new TBQSR_cellQuality[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(initialError, pmdObject);
	}

	//Do we also consider the effect of the positon in read (cycle)?
	considerPosition = params.parameterExists("considerPosition");
	if(considerPosition){
		maxPos = params.getParameterInt("maxPos");
		if(maxPos < 1) throw "Max position has to be larger than zero!";
		BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];
		for(int i=0; i<numReadGroups; ++i){
			BQSR_cells_readGroup_position[i] = new TBQSR_cellPosition[maxPos];
			for(int p=0; p<maxPos; ++p) BQSR_cells_readGroup_position[i][p].init(BQSR_cells_readGroup_quality, qualityIndex, pmdObject);
		}
	} else BQSR_cells_readGroup_position = NULL;

	//Do we also consider the context (dinucleotide)?
	considerContext = params.parameterExists("considerContext");
	if(considerContext){
		BQSR_cells_quality_context = new TBQSR_cellContext*[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q){
			BQSR_cells_quality_context[q] = new TBQSR_cellContext[numContexts]; //there are 25 contexts!
			for(int c=0; c<numContexts; ++c){
				if(considerPosition) BQSR_cells_quality_context[q][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex, pmdObject);
				else BQSR_cells_quality_context[q][c].init(BQSR_cells_readGroup_quality, qualityIndex, pmdObject);
			}
		}
	} else BQSR_cells_quality_context = NULL;
}

void TRecalibrationBQSR::addSite(TSite & site){
	if(site.referenceBase != 'N'){
		Base refBase = site.getRefBaseAsEnum();
		if(!qualityConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				BQSR_cells_readGroup_quality[(*it)->readGroup][qualityIndex->getIndex((*it)->quality)].addBase(*it, refBase);
			}
		}
		else if(considerPosition && !positionConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInRead >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position[(*it)->readGroup][(*it)->posInRead].addBase(*it, refBase);
			}
		} else if(considerContext && !contextConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				BQSR_cells_quality_context[qualityIndex->getIndex((*it)->quality)][(*it)->context].addBase(*it, refBase);
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(){
	//estimate epsilon for quality and readgroup, if not yet done
	estimationConverged = true;
	int numCellsNotConverged;
	int totCells;
	double maxF = 0.0;

	//readGroup x quality
	if(!qualityConverged){
		logfile->listFlush("Estimating epsilon for readGroup x quality table ...");
		qualityConverged = true;
		numCellsNotConverged = 0; totCells = 0;
		for(int i=0; i<numReadGroups; ++i){
			for(int q=0; q<qualityIndex->numQ; ++q){
				++totCells;
				if(!BQSR_cells_readGroup_quality[i][q].estimate(convergenceThreshold)){
					qualityConverged = false;
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_quality[i][q].F > maxF) maxF = BQSR_cells_readGroup_quality[i][q].F;
				}
			}
		}
		logfile->write(" done!");
		if(!qualityConverged){
			//empty all cells of BQSR table
			for(int i=0; i<numReadGroups; ++i){
				for(int q=0; q<qualityIndex->numQ; ++q){
					BQSR_cells_readGroup_quality[i][q].empty();
				}
			}

			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * qualityIndex->numQ));

			percent = 100 * (double) numCellsNotConverged / (double)totCells;

			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
			logfile->conclude("Largest F = " + toString(maxF));

			//set status
			estimationConverged = false;
			return estimationConverged;
		} else {
			logfile->list("Estimation converged in all cells!");
		}
	}

	//estimate epsilon for position, if not yet done
	if(considerPosition && !positionConverged){
		logfile->listFlush("Estimating epsilon for readGroup x position table ...");
		positionConverged = true;
		numCellsNotConverged = 0;
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(convergenceThreshold)){
					positionConverged = false;
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_position[i][p].F > maxF) maxF = BQSR_cells_readGroup_position[i][p].F;
				}
			}
		}
		if(!positionConverged){
			//empty all cells of BQSR table
			for(int i=0; i<numReadGroups; ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position[i][p].empty();
				}
			}
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
			logfile->conclude("Largest F = " + toString(maxF));

			//set status
			estimationConverged = false;
			return estimationConverged;
		} else {
			logfile->list("Estimation converged in all cells!");
		}
	}

	//estimate epsilon for context
	if(considerPosition && !positionConverged){
		logfile->listFlush("Estimating epsilon for quality x context table ...");
		contextConverged = true;
		numCellsNotConverged = 0;
		for(int q=0; q<qualityIndex->numQ; ++q){
			for(int c=0; c<numContexts; ++c){
				if(!BQSR_cells_quality_context[q][c].estimate(convergenceThreshold)){
					contextConverged = false;
					++numCellsNotConverged;
					if(BQSR_cells_quality_context[q][c].F > maxF) maxF = BQSR_cells_quality_context[q][c].F;
				}
			}
		}
		if(!contextConverged){
			//empty all cells of BQSR table
			for(int q=0; q<qualityIndex->numQ; ++q){
				for(int c=0; c<numContexts; ++c){
					BQSR_cells_quality_context[q][c].empty();
				}
			}

			int percent = 100.0 * ((double) numCellsNotConverged / (double) (qualityIndex->numQ * numContexts));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
			logfile->conclude("Largest F = " + toString(maxF));

			//set status
			estimationConverged = false;
			return estimationConverged;
		} else {
			logfile->list("Estimation converged in all cells!");
		}
	}

	//return true on final convergence
	return estimationConverged;
}

void TRecalibrationBQSR::writeToFile(std::string filenameTag){
	//write readGroup x Quality table
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	std::ofstream out(filename.c_str());
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int q=0; q<qualityIndex->maxQ + 1; ++q){
			out << it->ID << "\t" << q << "\tM\t" << BQSR_cells_readGroup_quality[i][qualityIndex->getIndex(q)].curEstimate << "\t" << BQSR_cells_readGroup_quality[i][q].numObservations << "\n";
		}
	}
	out.close();

	//write readGroup x position table
	if(considerPosition){
		filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
		std::ofstream out(filename.c_str());
		out << "ReadGroup\tPosition\tEventType\tEmpiricalQuality\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int i=0; i<numReadGroups; ++i, ++it){
			for(int p=0; p<maxPos; ++p){
				out << it->ID << "\t" << p << "\tM\t" << BQSR_cells_readGroup_position[i][p].curEstimate << "\t" << BQSR_cells_readGroup_position[i][p].numObservations << "\n";
			}
		}
		out.close();
	}

	//write readGroup x position table
	if(considerContext){
		filename = filenameTag + "_BQSR_Quality_Context_Table.txt";
		std::ofstream out(filename.c_str());
		out << "ReadGroup\tPosition\tEventType\tEmpiricalQuality\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int q=0; q<qualityIndex->maxQ + 1; ++q){
			for(int c=0; c<numContexts; ++c){
				out << it->ID << "\t" << q << "\tM\t" << BQSR_cells_quality_context[q][c].curEstimate << "\t" << BQSR_cells_quality_context[q][c].numObservations << "\n";
			}
		}
		out.close();
	}
}


