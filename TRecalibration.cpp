/*
 * TRecalibration.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: wegmannd
 */

#include "TRecalibration.h"

//---------------------------------------------------------------
//RecalibrationEMSite
//---------------------------------------------------------------
TRecalibrationEMSite::TRecalibrationEMSite(){

	std::cout << "CONST SIMPLE -----------------" << std::endl;

	initialized = false;
	q = NULL;
	context = NULL;
	readGroup = NULL;
	readGroupShifts = NULL;
	D = NULL;
	B = NULL;
	epsilon = NULL;
	P_g_given_d_oldBeta = NULL;
	numReads = 0;
};

TRecalibrationEMSite::TRecalibrationEMSite(TSite & site){
	numReads = site.bases.size();
	q = new double*[numReads];
	D = new double*[4];
	B = new double*[4];
	for(int g=0; g<4; ++g){
		D[g] = new double[numReads];
		B[g] = new double[numReads];
	}

	context = new int[numReads];
	epsilon = new double[numReads];
	readGroup = new int[numReads];
	readGroupShifts = new int[numReads];
	P_g_given_d_oldBeta = new double[4];
	initialized = true;
	int k=0;
	double epsilon;
	for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it, ++k){
		readGroup[k] = (*it)->readGroup;
		readGroupShifts[k] = readGroup[k] * 25; //shift by num params per read group!
		q[k] = new double[4];

		//we will work with the following q_ikl:
		// - transformed quality
		// - square of transformed quality
		epsilon = dePhred((*it)->quality);
		q[k][0] = log(epsilon / (1.0 - epsilon));
		q[k][1] = q[k][0] * q[k][0];

		// - position
		// - square of position
		q[k][2] = (*it)->posInRead;
		q[k][3] = q[k][2] * q[k][2];

		// - 20 context indicators (either 0.0 or 1.0)
		//only store which is one!
		context[k] = (*it)->context;

		//now also store D: D[ref][read]
		switch((*it)->getBaseAsEnum()){
			case A: D[0][k] = 0.0; //geno = AA
					D[1][k] = 1.0; //geno = CC
					D[2][k] = 1.0 - (*it)->PMD_GA; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case C: D[0][k] = 1.0; //geno = AA
					D[1][k] = (*it)->PMD_CT; //geno = CC
					D[2][k] = 1.0; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case G: D[0][k] = 1.0; //geno = AA
					D[1][k] = 1.0; //geno = CC
					D[2][k] = (*it)->PMD_GA; //geno = GG
					D[3][k] = 1.0; //geno = TT
					break;
			case T: D[0][k] = 1.0; //geno = AA
					D[1][k] = 1.0 - (*it)->PMD_CT; //geno = CC
					D[2][k] = 1.0; //geno = GG
					D[3][k] = 0.0; //geno = TT
					break;
			case N:
					D[0][k] = 0.0;
					D[1][k] = 0.0;
					D[2][k] = 0.0;
					D[3][k] = 0.0;
					break;
		}

		//now store B
		for(int g=0; g<4; ++g){
			B[g][k] = 4.0 / 3.0 * D[g][k] - 1.0;
		}
	}
};

TRecalibrationEMSite::~TRecalibrationEMSite(){
	if(initialized){
		for(int i=0; i<4; ++i){
			delete[] D[i];
			delete[] B[i];
		}
		for(int i=0; i<numReads; ++i){
			delete[] q[i];
		}
		delete[] q;
		delete[] D;
		delete[] B;
		delete[] context;
		delete[] readGroup;
		delete[] readGroupShifts;
		delete[] epsilon;
		delete[] P_g_given_d_oldBeta;
	}
}

void TRecalibrationEMSite::calcEpsilon(double** params){
	//calc epsilon using parameter estimates provided
	double eta;
	for(int k=0; k<numReads; ++k){
		eta = params[readGroup[k]][0];
		for(int p=0; p<4; ++p){ //loop over all parameters except beta0
			eta += params[readGroup[k]][p+1] * q[k][p];
		}
		//eta += params[readGroup[k]][context[k]];

		if(eta > 22.2) epsilon[k] = 0.9999999999;
		else if(eta < -23.02685) epsilon[k] = 0.0000000001;
		else {
			eta = exp(eta);
			epsilon[k] = eta / (1.0 + eta);
		}
	}
}

double TRecalibrationEMSite::fill_P_g_given_d_beta_AND_calcLL(double** oldParams, double* freqs){
	calcEpsilon(oldParams);

	//over all genotypes
	double P_g_given_d_theta_denominator = 0.0;
	double tmp;
	for(int g=0; g<4; ++g){
		tmp = 1.0;
		//loop over all reads
		for(int k=0; k<numReads; ++k){
			tmp *= B[g][k] * epsilon[k] - D[g][k] + 1;
		}
		P_g_given_d_oldBeta[g] = tmp * freqs[g];
		P_g_given_d_theta_denominator += P_g_given_d_oldBeta[g];

	}

	//calculate P(g|d, theta)
	for(int g=0; g<4; ++g){
		P_g_given_d_oldBeta[g] = P_g_given_d_oldBeta[g] / P_g_given_d_theta_denominator;
	}

	//return LL = P_g_given_d_theta_denominator
	return log(P_g_given_d_theta_denominator);
}

double TRecalibrationEMSite::calcQ(double** newParams){
	calcEpsilon(newParams);

	//now calculate P(d, g, new params)
	double P_d_given_g_beta;
	double Q = 0.0;

	for(int g=0; g<4; ++g){
		P_d_given_g_beta = 1.0;
		//loop over all reads
		for(int k=0; k<numReads; ++k){
			P_d_given_g_beta *= B[g][k] * epsilon[k] - D[g][k] + 1;
		}

		if(P_d_given_g_beta < 1.0E-50) P_d_given_g_beta = 1.0E-50;
		Q += P_g_given_d_oldBeta[g] * log(P_d_given_g_beta);
	}

	return Q;
}

void TRecalibrationEMSite::addToJacobianAndF(arma::mat & Jacobian, arma::vec & F, double** params){
	//calculate epsilon with current parameters
	calcEpsilon(params);

	//tmp variables
	double* weights = new double[numReads];
	double* eps1MinusEps = new double[numReads];
	double* oneMinus2Eps = new double[numReads];
	double* weightJacobian = new double[numReads];
	for(int k=0; k<numReads; ++k){
		eps1MinusEps[k] = epsilon[k] * (1.0 - epsilon[k]);
		oneMinus2Eps[k] = 1.0 - 2.0 * epsilon[k];
	}
	double tmp;
	int tmpIndex;

	//DEBUG-------------------------------------------------------------------
	//set all loops to go over the first two params only
	//DEBUG-------------------------------------------------------------------

	//fill F and Jacobian
	for(int g=0; g<4; ++g){
		//calc weights
		//------------

		for(int k=0; k<numReads; ++k){
			weights[k] = B[g][k] / (1.0 - D[g][k] + B[g][k] * epsilon[k]) * eps1MinusEps[k];
			weightJacobian[k] = P_g_given_d_oldBeta[g] * weights[k] * (oneMinus2Eps[k] - weights[k]);
		}

		//add to F
		//--------
		//beta 0
		for(int k=0; k<numReads; ++k){
			F(readGroupShifts[k]) += P_g_given_d_oldBeta[g] * weights[k];
		}

		//others
		for(int k=0; k<numReads; ++k){
			tmp = P_g_given_d_oldBeta[g] * weights[k];
			//now all 4 covariates except context. Derivatives are given by the q's
			tmpIndex = 1 + readGroupShifts[k];
			for(int m=0; m<4; ++m){ //loop over all parameters except beta0
				F(m + tmpIndex) += tmp * q[k][m];
			}
			//now context: start at position 5 in F!
			//F(context[k] + 5 + readGroupShifts[k]) += tmp;
		}

		//add to Jacobian
		//---------------
		for(int k=0; k<numReads; ++k){
			tmp = weightJacobian[k];

			//beta0
			Jacobian(readGroupShifts[k],readGroupShifts[k]) += tmp;

			//first row
			for(int m=0; m<4; ++m){ //loop over all parameters except beta0 and context
				Jacobian(readGroupShifts[k],m + 1 + readGroupShifts[k]) += tmp * q[k][m];
			}

			//all other rows except context
			tmpIndex = readGroupShifts[k] + 1;
			for(int row=0; row<4; ++row){ //loop over all parameters except beta0
				for(int col=row; col<4; ++col){ //loop over all parameters except beta0
					Jacobian(tmpIndex + row, tmpIndex + col) +=  tmp * q[k][col] * q[k][row];
				}
			}


			/*
			//context column: first five rows
			tmpIndex = readGroupShifts[k] + context[k] + 5;
			Jacobian(readGroupShifts[k], tmpIndex) += tmp;
			for(int p=0; p<4; ++p){
				Jacobian(readGroupShifts[k] + p + 1, tmpIndex) += tmp * q[k][p];
			}

			//context x context: only add to diagonal, as all others are 0
			Jacobian(tmpIndex, tmpIndex) += tmp;
			*/
		}
	} //end loop over genotypes

	//delete tmp variables
	delete[] weights;
	delete[] eps1MinusEps;
	delete[] oneMinus2Eps;
	delete[] weightJacobian;
}

//---------------------------------------------------------------
//TRecalibrationEMWindow
//---------------------------------------------------------------
TRecalibrationEMWindow::TRecalibrationEMWindow(TBaseFrequencies* baseFreqs){
	freqs = new double[4];
	for(int i=0; i<4; ++i) freqs[i] = (*baseFreqs)[i];
}

void TRecalibrationEMWindow::addSite(TSite & site){
	sites.push_back(new TRecalibrationEMSite(site));
}

double TRecalibrationEMWindow::fill_P_g_given_d_beta_AND_calcLL(double** oldParams){
	double LL = 0.0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		LL += (*site)->fill_P_g_given_d_beta_AND_calcLL(oldParams, freqs);
	}
	return LL;
}

double TRecalibrationEMWindow::calcQ(double** newParams){
	double Q = 0.0;
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		Q += (*site)->calcQ(newParams);
	}
	return Q;
}

void TRecalibrationEMWindow::addToJacobianAndF(arma::mat & Jacobian, arma::vec & F, double** params){
	for(std::vector<TRecalibrationEMSite*>::iterator site = sites.begin(); site != sites.end(); ++site){
		(*site)->addToJacobianAndF(Jacobian, F, params);
	}
}
//---------------------------------------------------------------
//TRecalibrationEM
//---------------------------------------------------------------
TRecalibrationEM::TRecalibrationEM(BamTools::SamHeader* BamHeader, TParameters & args, TLog* Logfile){
	//create data structure to store q_ikl for each observation
	//we will work with the following q_ikl (per read group):
	// - transformed quality
	// - square of transformed quality
	// - position
	// - square of position
	// - 20 context indicators (either 0.0 or 1.0)
	// -> in total, 25 variables to estimate (incl. intercept at last position)
	//if these are changed, TRecalibrationEMSite needs to be changed!
	numParams = 5;

	//rad groups and log file
	bamHeader = BamHeader;
	numReadGroups = bamHeader->ReadGroups.Size();
	readGroupNames = new std::string[numReadGroups];
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int r=0; r<numReadGroups; ++r, ++it){
		readGroupNames[r] = it->ID;
	}
	totNumParams = numParams * numReadGroups;

	logfile = Logfile;

	//create parameter arrays
	params = new double*[numReadGroups];
	newParams = new double*[numReadGroups];
	tmpParams = new double*[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		params[r] = new double[numParams];
		newParams[r] = new double[numParams];
		tmpParams[r] = new double[numParams];
	}

	//Are the values provided?
	estimatetionRequired = false;
	if(args.parameterExists("recal")){
		//read parameters from file

		//TODO!


		numSitesAdded = 0; //SET TO WHAT IT WILL BE!

		//check if we anyway estimate things
		if(args.parameterExists("estimateRecal")) estimatetionRequired = true;
	} else {
		estimatetionRequired = true;
		//set initial values: all to 0 except beta0 (quality) = 1
		for(int r=0; r<numReadGroups; ++r){
			for(int i=0; i<numParams; ++i){
				params[r][i] = 0.0;
				newParams[r][i] = 0.0;
			}

			//debug: no need to set [0]!
			params[r][0] = 0.0;
			newParams[r][0] = 0.0;

			params[r][1] = 1.0;
			newParams[r][1] = 1.0;

		}
		numSitesAdded = 0;
	}

	//read estimation parameters, if required
	if(estimatetionRequired){
		logfile->startIndent("Will run EM to estimate recalibration parameters:");
		numEMIterations = args.getParameterIntWithDefault("iterations", 10);
		logfile->list("Will perform at max " + toString(numEMIterations) + " EM iterations.");
		maxEpsilon = args.getParameterDoubleWithDefault("maxEps", 0.000001);
		logfile->list("Will stop EM when deltaLL < " + toString(maxEpsilon));
		NewtonRaphsonNumIterations = args.getParameterIntWithDefault("NRiterations", 10);
		logfile->list("Will conduct at max " + toString(NewtonRaphsonNumIterations) + " Newton-Raphson iterations");
		NewtonRaphsonMaxF = args.getParameterDoubleWithDefault("maxF", 0.0001);
		logfile->list("Will stop Newton-Raphson when F < " + toString(NewtonRaphsonMaxF));
		logfile->endIndent();

		//initialize vriables for EM
		Jacobian.resize(totNumParams, totNumParams);
		Jacobian.zeros();
		F.resize(totNumParams);
		F.zeros();
		JxF.resize(totNumParams, totNumParams);
		JxF.zeros();
	} else {
		numEMIterations = -1;
		maxEpsilon = 0.0;
		NewtonRaphsonNumIterations = -1;
		NewtonRaphsonMaxF = 0.0;
	}
}


void TRecalibrationEM::addNewWindow(TBaseFrequencies* freqs){
	windows.push_back(new TRecalibrationEMWindow(freqs));
	//windows.emplace_back(freqs);
	//set iterator
	curWindow = windows.end(); --curWindow;
}

void TRecalibrationEM::addSite(TSite & site){
	(*curWindow)->addSite(site);
	++numSitesAdded;
}

double TRecalibrationEM::getErrorRate(TBase* base, double** theseParams){
	//eta = beta0 + SUM_i beta[i] * q[i]
	int rg = base->readGroup;
	double eta = theseParams[rg][0];

	// q[1] is transformed quality
	double tmp = dePhred(base->quality);
	tmp = log(tmp / (1.0 + tmp));
	eta += theseParams[rg][1] * tmp;

	//q[2] is square of transformed quality
	eta += theseParams[rg][2] * tmp * tmp;

	//q[3] is position
	eta += theseParams[rg][3] * base->posInRead;

	//q[4] is square of position
	eta += theseParams[rg][4] * base->posInRead * base->posInRead;

	//q[5] until q[24] are indicators for the context. Just pick the matching one!
	//eta += theseParams[rg][base->context];

	//now calculate epsilon from eta
	if(eta > 22.2) return 0.9999999999;
	if(eta < -23.02685) return 0.0000000001;

	eta = exp(eta);
	return eta / (1.0 + eta);
}

double TRecalibrationEM::getErrorRate(TBase* base){
	return getErrorRate(base, params);
}

void TRecalibrationEM::runNewtonRaphson(double** theseParams, int & maxNewtonRaphsonIteratios, double & maxFThreshold, TLog* logfile, std::string debugFilename){
	//variables
	double maxF;
	int index;
	double lambda; //used in backtracking
	bool acceptMove;
	bool NRconverged = false;

	//calculate Q at current location
	double Q;
	double curQ = 0.0;
	for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
		curQ += (*curWindow)->calcQ(theseParams);
	}

	//open debug file
	std::ofstream out(debugFilename.c_str());
	if(!out) throw "Failed to open output file '" + debugFilename + "'!";
	//add header
	out << "iteration";
	for(int i=0; i<numParams; ++i) out << "\tbeta'" << i;
	for(int i=0; i<numParams; ++i) out << "\tF" << i;
	for(int i=0; i<numParams; ++i) out << "\tbeta" << i;
	out << std::endl;

	//run up to maxNewtonRaphsonIteratios iterations, but stop if max(F) < maxFThreshold
	logfile->startIndent("Running Newton-Raphson optimization:");
	for(int i=0; i<maxNewtonRaphsonIteratios; ++i){
		logfile->startIndent("Running iteration " + toString(i+1) + ":");
		logfile->listFlush("Calculating Jacobian and graident ...");
		out << i;

		//set to zero
		Jacobian.zeros();
		F.zeros();

		//fill Jacobin and F: loop over all sites
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			(*curWindow)->addToJacobianAndF(Jacobian, F, theseParams);
		}

		//Need to copy numbers to other triangle in Jacobian, as only upper triangle is filled when parsing sites
		for(int i=0; i<(totNumParams-1); ++i){
			for(int j=i+1; j<totNumParams; ++j){
				//copy from upper triangle to lower triangle
				Jacobian(j,i) = Jacobian(i,j);
			}
		}

		//scale F and J by 1/#sites
		Jacobian = Jacobian / (double) numSitesAdded;
		F = F / (double) numSitesAdded;

		//now calculate J^-1 x F

		/*
		std::cout << std::endl << "-------JACOBIAN-------------------------------" << std::endl;
		std::cout << Jacobian << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "F: " << F << std::endl;
		std::cout << "----------------------------------------------" << std::endl;
		std::cout << "det(J) = " << det(Jacobian) << std::endl;
		std::cout << std::endl;
		std::cout << "----------------------------------------------" << std::endl;

*/
		//print beta'
		for(int i=0; i<numParams; ++i){
			out << "\t" << theseParams[0][i];
		}

		//print F
		for(int i=0; i<numParams; ++i){
			out << "\t" << F[i];
		}


		if(solve(JxF, Jacobian, F)){
			logfile->write(" done!");

			/*
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << "JxF " << JxF << std::endl;
			std::cout << "----------------------------------------------" << std::endl;
*/
			//update params for each read group using backtracking
			lambda = 1.0;
			acceptMove = false;
			while(!acceptMove){
				logfile->listFlush("Proposing move with lambda = " + toString(lambda) + " ...");
				//estimate new params
				for(int r=0; r<numReadGroups; ++r){
					index = r*numParams;
					for(int i=0; i<numParams; ++i){
						tmpParams[r][i] = theseParams[r][i] - lambda * JxF(index + i);
					}
				}

				//calculate Q at new location
				Q = 0.0;
				for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
					Q += (*curWindow)->calcQ(tmpParams);
				}

				//check if we accept or backtrack
				if(Q > curQ){
					acceptMove = true; //accept
					logfile->write(" accepting move!");
					logfile->conclude("Q was reduced from " + toString(curQ) + " to " + toString(Q));
					curQ = Q;
					//store new params
					for(int r=0; r<numReadGroups; ++r){
						for(int i=0; i<numParams; ++i){
							theseParams[r][i] = tmpParams[r][i];
						}
					}
				}
				else{
					lambda = lambda / 2.0; //backtrack;
					logfile->write(" rejecting move!");
					if(lambda < 0.000000001){
						acceptMove = true; //accept
						NRconverged = true;
						logfile->conclude("No improvement even with lambda = " + toString(lambda) + ", aborting Newton-Raphson.");
					}
				}
			}
		} else throw "Issue solving JxF in TRecalibrationEM::runNewtonRalphson()!";

		//print beta
		for(int i=0; i<numParams; ++i){
			out << "\t" << theseParams[0][i];
		}
		out << std::endl;

		//get largest gradient (F) to check if we break
		maxF = 0.0;
		for(int i=0; i<numParams; ++i){
			if(fabs(F(i)) > maxF) maxF = fabs(F(i));
		}
		logfile->conclude("max(F) = " + toString(maxF));
		logfile->endIndent();
		if(maxF < maxFThreshold || NRconverged) break;
	}
	out.close();
	logfile->endIndent();
}

void TRecalibrationEM::runEM(std::string outputName){
	logfile->startNumbering("Running EM algorithm to find MLE recalibration parameters:");
	double LL, deltaLL, oldLL = 0.0;
	std::ofstream out;
	std::string filename;

	//running iterations
	for(int iter = 0; iter < numEMIterations; ++iter){
		logfile->number("EM Iteration:"); logfile->addIndent();

		//calculate P(g|d, oldbeta) for all sites and calculate LL
		LL = 0.0;
		logfile->listFlush("Calculating P(g|d, beta') ...");
		for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
			LL += (*curWindow)->fill_P_g_given_d_beta_AND_calcLL(params);
		}
		logfile->write(" done!");
		logfile->conclude("Current Log Likelihood = " + toString(LL));
		deltaLL = LL - oldLL;
		logfile->conclude("Epsilon = " + toString(deltaLL));


		//DEBUG--------------------------------------------------------
		//calc Q surface for current old params
		//calcQSurface(outputName + "_Qsurface_EMiteration_" + toString(iter) + ".txt", 21);
		//DEBUG--------------------------------------------------------

		//fill vector of new params by copying current values
		for(int r=0; r<numReadGroups; ++r){
			for(int i=0; i<numParams; ++i){
				newParams[r][i] = params[r][i];
			}
		}

		//check if we break based on LL
		if(iter > 0 && deltaLL < maxEpsilon){
			logfile->conclude("EM has converged (epsilon < " + toString(maxEpsilon) + ")");
			break;
		}
		else oldLL = LL;

		//run NewtonRaphson until convergence
		runNewtonRaphson(newParams, NewtonRaphsonNumIterations, NewtonRaphsonMaxF, logfile, outputName + "_NewtonRaphson_" + toString(iter) + ".txt");

		//save parameters
		for(int r=0; r<numReadGroups; ++r){
			for(int i=0; i<numParams; ++i){
				params[r][i] = newParams[r][i];
			}
		}

		//write current estimates to file
		filename = outputName + "_recalibrationEM_Loop" + toString(iter) + ".txt";
		logfile->listFlush("Writing current estimates to file '" + filename + "' ...");
		writeCurrentEstimates(filename, LL);
		logfile->write(" done!");

		//end loop
		logfile->endIndent();
	}

	//finalize
	logfile->endNumbering();

	//writing final estimates
	filename = outputName + "_recalibrationEM.txt";
	logfile->list("Writing final estimates to file '" + filename + "' ...");
	writeCurrentEstimates(filename, LL);
	logfile->write(" done!");

	//calc LL surface
	//calcLikelihoodSurface(outputName + "_LLsurface.txt", 21);
}

void TRecalibrationEM::writeCurrentEstimates(std::string filename, double & LL){
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	writeHeader(out);
	writeParams(out, LL);
	out.close();
}

void TRecalibrationEM::writeHeader(std::ofstream & out){
	out << "readGroup";
	for(int i=0; i<numParams; ++i)
		out << "\tbeta" << i;
	out << "\tLL" << std::endl;
}

void TRecalibrationEM::writeParams(std::ofstream & out, double & LL){
	for(int r=0; r<numReadGroups; ++r){
		out << readGroupNames[r];
		for(int i=0; i<numParams; ++i){
			out << "\t" << params[r][i];
		}
		out << "\t" << LL;
		out << std::endl;
	}
}

void TRecalibrationEM::calcLikelihoodSurface(std::string filename, int numMarginalGridPoints){
	double LL;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tLL" << std::endl;

	//set min, max and step for each parameter
	double min[5];
	min[0] = -5.0;
	min[1] = -5.0;
	min[2] = -1.0;
	min[3] = -1.0;
	min[4] = -1.0;


	double max[5];
	max[0] = 10.0;
	max[1] = 10.0;
	max[2] = 1.0;
	max[3] = 1.0;
	max[4] = 1.0;

	double step[5];
	for(int i=0; i<5; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//without last two
	for(int r=0; r<numReadGroups; ++r){
		params[r][3] = 0.0;
		params[r][4] = 0.0;
	}

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		for(int r=0; r<numReadGroups; ++r) params[r][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			for(int r=0; r<numReadGroups; ++r) params[r][1] = min[1] + p2 * step[1];
			for(int p3=0; p3<numMarginalGridPoints; ++p3){
				for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate LL
						LL = 0.0;
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							LL += (*curWindow)->fill_P_g_given_d_beta_AND_calcLL(params);
						}

						//write to file
						for(int i=0; i<5; ++i) out << params[0][i] << "\t";
						out << LL << std::endl;
					//}
				//}
			}
		}
	}

	//close file
	out.close();
}


void TRecalibrationEM::calcQSurface(std::string filename, int numMarginalGridPoints){
	double Q;

	//open outputfile
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open output file '" + filename + "'!";
	out << "beta0\tbeta1\tQ" << std::endl;

	//set min, max and step for each parameter
	double min[2];
	min[0] = -5.0;
	min[1] = -5.0;

	double max[2];
	max[0] = 10.0;
	max[1] = 10.0;

	double step[2];
	for(int i=0; i<2; ++i){
		step[i] = (max[i] - min[i]) / (numMarginalGridPoints - 1.0);
	}

	//print old params

	//Loop over parameters
	for(int p1=0; p1<numMarginalGridPoints; ++p1){
		for(int r=0; r<numReadGroups; ++r) newParams[r][0] = min[0] + p1 * step[0];
		for(int p2=0; p2<numMarginalGridPoints; ++p2){
			for(int r=0; r<numReadGroups; ++r) newParams[r][1] = min[1] + p2 * step[1];
			//for(int p3=0; p3<numMarginalGridPoints; ++p3){
				//for(int r=0; r<numReadGroups; ++r) params[r][2] = min[2] + p3 * step[2];
				//for(int p4=0; p4<numMarginalGridPoints; ++p4){
					//for(int r=0; r<numReadGroups; ++r) params[r][3] = min[3] + p4 * step[3];
					//for(int p5=0; p5<numMarginalGridPoints; ++p5){
						//for(int r=0; r<numReadGroups; ++r) params[r][4] = min[4] + p5 * step[4];


						//calculate Q
						Q = 0.0;
						//logfile->listFlush("Calculating Q at {" + toString(params[0][0]) + ", " + toString(params[0][1]) + "} ...");
						for(curWindow = windows.begin(); curWindow != windows.end(); ++curWindow){
							Q += (*curWindow)->calcQ(newParams);
						}
						//logfile->write(" done!");
						//logfile->conclude("Current Q = " + toString(Q));

						//write to file
						for(int i=0; i<2; ++i) out << newParams[0][i] << "\t";
						out << Q << std::endl;
					//}
				//}
			//}
		}
	}

	//close file
	out.close();
}

//---------------------------------------------------------------
//TRecalibrationBQSR_cell BQSR
//---------------------------------------------------------------
TBQSR_cell::TBQSR_cell(){
	curEstimate = 0.0;
	estimationConverged = false;
	firstDerivative = 0.0;
	firstDerivativeSave = 0.0;
	secondDerivative = 0.0;
	secondDerivativeSave = 0.0;
	numObservations = 0.0;
	numObservationsTmp = 0.0;
	numMatches = 0.0;
	F = 0.0;
	LL = 0.0;
}

void TBQSR_cell::init(double initialError){
	curEstimate = initialError;
}

void TBQSR_cell::set(double error, std::string & NumObservations){
	curEstimate = error;
	if(NumObservations == "-") numObservations = 0;
	else numObservations = pow(10.0, stringToDouble(NumObservations));
};

void TBQSR_cell::empty(){
	if(!estimationConverged){
		numObservationsTmp = 0;
		numMatches = 0;
		firstDerivativeSave = firstDerivative;
		secondDerivativeSave = secondDerivative;
		firstDerivative = 0.0;
		secondDerivative = 0.0;
		LL = 0.0;
	}
}

void TBQSR_cell::reopenEstimation(){
	estimationConverged = false;
	empty();
}

double TBQSR_cell::getD(TBase* base, Base & RefBase){
	double D = 0.0;
	switch(base->getBaseAsEnum()){
		case A: if(RefBase == A){
					D = 1.0;
					break;
				}
				if(RefBase == G) D = base->PMD_GA;
				break;
		case C: if(RefBase == C) D = 1.0 - base->PMD_CT;
				break;
		case G: if(RefBase == G) D = 1.0 - base->PMD_GA;
				break;
		case T: if(RefBase == C) D = base->PMD_CT;
		        else if(RefBase == T) D = 1.0;
				break;
		case N: throw "Can not add base with unknown reference to BQSR cell!";
	}
	return D;
}

void TBQSR_cell::addToDerivatives(TBase* base, Base & RefBase){
	double D = getD(base, RefBase);

	double oneMinus4D = 1.0 - 4.0 * D;
	firstDerivative += oneMinus4D / (-4.0*D*curEstimate + 3.0*D + curEstimate);
	double tmpF = oneMinus4D / (D*(3.0-4.0*curEstimate) + curEstimate);
	secondDerivative -= tmpF * tmpF;

	++numObservationsTmp;
}


void TBQSR_cell::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		addToDerivatives(base, RefBase);
		if(base->getBaseAsEnum() == RefBase) ++numMatches;
	}
}

void TBQSR_cell::runNewtonRaphson(double & convergenceThreshold){
	curEstimate = curEstimate - firstDerivative / secondDerivative;
	//decide on convergence
	F = fabs(firstDerivative / numObservations);
	if(F < convergenceThreshold) estimationConverged = true;
}

bool TBQSR_cell::estimate(double & convergenceThreshold, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
		} else if(numMatches >= numObservations){ //epsilon = 0
			curEstimate = 0.0;
			estimationConverged = true;
		} else if(numMatches < 1.0){ // epsilon = 1.0
			curEstimate = 1.0;
			estimationConverged = true;
		} else {
			//need Newton-Raphson to estimate epsilon
			double oldEstimate = curEstimate;
			runNewtonRaphson(convergenceThreshold);

			//check boundaries
			if(curEstimate < 0.0){
				curEstimate = 0.000000001;
				if(oldEstimate == 0.00000001)
					estimationConverged = true; //if estimate is repeatedly below, accept
			} else if(curEstimate > 1.0){
				curEstimate = 0.999999999;
				if(oldEstimate == 0.999999999)
					estimationConverged = true; //if estimate is repeatedly above, accept
			}

			//do not allow big jump in quality -> max +/- 10!
			if(curEstimate / oldEstimate > 10.0){
				curEstimate = oldEstimate * 10.0;
			} else if(oldEstimate / curEstimate > 10.0){
				curEstimate = oldEstimate / 10.0;
			}
		}
	}
	return estimationConverged;
}

std::string TBQSR_cell::getNumObsForPrinting(){
	if(numObservations == 0) return "-";
	else return toString(log10(numObservations));
}
//---------------------------------------------------------------
TBQSR_cellPosition::TBQSR_cellPosition():TBQSR_cell(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
}

void TBQSR_cellPosition::init(TBQSR_cell** gotBQSR_cells_quality_readGroup, TQualityIndex* QualityIndex){
	BQSR_cells_readGroup_quality = gotBQSR_cells_quality_readGroup;
	qualityIndex = QualityIndex;
}

void TBQSR_cellPosition::addToDerivatives(TBase* base, Base & RefBase, double & epsilon){
	double D = getD(base, RefBase);

	double epsMinus4Deps = epsilon - 4.0 * D * epsilon;
	firstDerivative += epsMinus4Deps / (-4.0*D*epsilon*curEstimate + 3.0*D + epsilon*curEstimate);
	double tmpF = epsMinus4Deps / (D*(3.0-4.0*epsilon*curEstimate) + epsilon*curEstimate);
	secondDerivative -= tmpF * tmpF;

	++numObservationsTmp;

	//LL += log((1.0-D) * epsilon*curEstimate / 3.0 + D * (1.0 - epsilon*curEstimate));
}

void TBQSR_cellPosition::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		addToDerivatives(base, RefBase, BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate);
	}
}

bool TBQSR_cellPosition::estimate(double & convergenceThreshold, long & minObservations){
	if(!estimationConverged){
		//set the number of observations this estimate was based on
		numObservations = numObservationsTmp;

		if(numObservations < minObservations){ //keep current estimate
			estimationConverged = true;
			return estimationConverged;
		} else {
			//need Newton-Raphson to estimate epsilon
			double oldEstimate = curEstimate;
			runNewtonRaphson(convergenceThreshold);

			//check boundaries
			if(curEstimate < 0.0){
				curEstimate = 0.01;
				if(oldEstimate == 0.01)
					estimationConverged = true; //if estimate is repeatedly below, accept
			} else if(curEstimate > 10000.0){
				curEstimate = 100.0;
				if(oldEstimate == 100.0)
					estimationConverged = true; //if estimate is repeatedly above, accept
			}
		}
	}
	return estimationConverged;
}

//---------------------------------------------------------------
TBQSR_cellPositionRev::TBQSR_cellPositionRev():TBQSR_cellPosition(){
	BQSR_cells_readGroup_quality = NULL;
	qualityIndex = NULL;
	curEstimate = 1.0;
	BQSR_cells_readGroup_position = NULL;
	considerPosition = false;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position = gotBQSR_quality_position;
	considerPosition = true;
}

void TBQSR_cellPositionRev::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex){
	TBQSR_cellPosition::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position = NULL;
}

void TBQSR_cellPositionRev::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
		if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
TBQSR_cellContext::TBQSR_cellContext():TBQSR_cellPositionRev(){
	BQSR_cells_readGroup_position_rev = NULL;
	considerPositionRev = false;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPosition** gotBQSR_quality_position, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, gotBQSR_quality_position, QualityIndex);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TBQSR_cellPositionRev** gotBQSR_quality_position_rev, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position_rev = gotBQSR_quality_position_rev;
	considerPositionRev = true;
}


void TBQSR_cellContext::init(TBQSR_cell** gotBQSR_quality_readGroup, TQualityIndex* QualityIndex){
	TBQSR_cellPositionRev::init(gotBQSR_quality_readGroup, QualityIndex);
	BQSR_cells_readGroup_position_rev = NULL;
}

void TBQSR_cellContext::addBase(TBase* base, Base & RefBase){
	if(!estimationConverged){
		double epsilonAlpha = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
		if(considerPosition) epsilonAlpha *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
		if(considerPositionRev) epsilonAlpha *= BQSR_cells_readGroup_position_rev[base->readGroup][base->posInReadRev].curEstimate;
		addToDerivatives(base, RefBase, epsilonAlpha);
	}
}

//---------------------------------------------------------------
//Recalibration BQSR
//---------------------------------------------------------------
TRecalibrationBQSR::TRecalibrationBQSR(BamTools::SamHeader* BamHeader, TParameters & params, TLog* Logfile){
	logfile = Logfile;
	bamHeader = BamHeader;
	numReadGroups = bamHeader->ReadGroups.Size();
	estimatetionRequired = false;
	estimationConverged = false;
	numContexts = 20;
	qualityIndex = NULL;
	maxPos = 0;

	//check if BQSR table readGroup x Quality is given, or has to be estimated
	initializeBQSRReadGroupQualityTable(params);

	//Do we also consider the effect of the position in read (cycle)?
	initializeBQSRReadGroupPositionTable(params);
	initializeBQSRReadGroupPositionReverseTable(params);

	//Do we also consider the context (dinucleotide)?
	initializeBQSRReadGroupContextTable(params);

	//read Newton-Raphson arguments from user
	convergenceThreshold = params.getParameterDoubleWithDefault("maxF", 0.0000001);
	if(estimatetionRequired) logfile->list("Stopping Newton-Raphson if F < " + toString(convergenceThreshold));

	//get minimal number of observations to conduct estimation
	minObservations = params.getParameterLongWithDefault("minObservations", 1000);
}

int TRecalibrationBQSR::findReadGroupIndex(std::string & name){
	int i = 0;
	for(std::vector<BamTools::SamReadGroup>::iterator it = bamHeader->ReadGroups.Begin(); it !=  bamHeader->ReadGroups.End(); ++it, ++i){
		if(name == it->ID) return i;
	}
	return -1;
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTable(TParameters & params){
	if(params.parameterExists("BQSRQuality")) initializeBQSRReadGroupQualityTableFromFile(params);
	else {
		qualityConverged = false;
		estimatetionRequired = true;
		int minQ = params.getParameterIntWithDefault("minQ", 1);
		int maxQ = params.getParameterIntWithDefault("maxQ", 100);
		logfile->list("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
		qualityIndex = new TQualityIndex(minQ, maxQ);

		//initialize BQSR table
		BQSR_cells_readGroup_quality = new TBQSR_cell*[numReadGroups];
		for(int i=0; i<numReadGroups; ++i){
			BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
			for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(dePhred(qualityIndex->getQuality(q)));
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupQualityTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRQuality");
	logfile->listFlush("Constructing BQSR readGroup x quality table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x quality table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_quality = new TBQSR_cell*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	int minQ = 100;
	int maxQ = 0;
	int q;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get min and max quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			q = stringToInt(vec[1]);
			if(q > maxQ) maxQ = q;
			if(q < minQ) minQ = q;
		}
	}

	//initialize quality index
	qualityIndex = new TQualityIndex(minQ, maxQ);

	//create corresponding objects
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_quality[i] = new TBQSR_cell[qualityIndex->numQ];
		for(int q=0; q<qualityIndex->numQ; ++q) BQSR_cells_readGroup_quality[i][q].init(dePhred(qualityIndex->getQuality(q)));
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double quality;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				q = stringToInt(vec[1]);
				quality = stringToDouble(vec[3]);
				BQSR_cells_readGroup_quality[readGroup][qualityIndex->getIndex(q)].set(dePhred(quality), vec[4]);
			}
		}
	}

	//set that no estimation is not required, unless asked for
	if(params.parameterExists("estimateBQSRQuality")) qualityConverged = false;
	else qualityConverged = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering qualities between " + toString(minQ) + " and " + toString(maxQ));
}


void TRecalibrationBQSR::initializeBQSRReadGroupPositionTable(TParameters & params){
	if(params.parameterExists("BQSRPosition")) initializeBQSRReadGroupPositionTableFromFile(params);
	else {
		positionConverged = false;
		considerPosition = params.parameterExists("estimateBQSRPosition");
		if(considerPosition){
			estimatetionRequired = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions up to " + toString(maxPos));
			BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];
			for(int i=0; i<numReadGroups; ++i){
				BQSR_cells_readGroup_position[i] = new TBQSR_cellPosition[maxPos];
				for(int p=0; p<maxPos; ++p) BQSR_cells_readGroup_position[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
			}
		} else {
			BQSR_cells_readGroup_position = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPosition");
	logfile->listFlush("Constructing BQSR readGroup x position table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position = new TBQSR_cellPosition*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_position[i] = new TBQSR_cellPosition[maxPos];
		isListed[i] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			BQSR_cells_readGroup_position[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
			isListed[i][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPosition")) positionConverged = false;
	else positionConverged = true;
	considerPosition = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering positions up to " + toString(maxPos));
}


//the functions are almost identical to the other position -> put in class!
void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTable(TParameters & params){
	if(params.parameterExists("BQSRPositionReverse")) initializeBQSRReadGroupPositionReverseTableFromFile(params);
	else {
		positionReverseConverged = false;
		considerPositionReverse = params.parameterExists("estimateBQSRPositionReverse");
		if(considerPositionReverse){
			estimatetionRequired = true;
			maxPos = params.getParameterInt("maxPos");
			if(maxPos < 1) throw "Max position has to be larger than zero!";
			logfile->list("Considering positions reverse up to " + toString(maxPos));
			BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[numReadGroups];
			for(int i=0; i<numReadGroups; ++i){
				BQSR_cells_readGroup_position_reverse[i] = new TBQSR_cellPositionRev[maxPos];
				for(int p=0; p<maxPos; ++p){
					if(considerPosition) BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
					else BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
				}
			}
		} else {
			BQSR_cells_readGroup_position_reverse = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupPositionReverseTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRPositionReverse");
	logfile->listFlush("Constructing BQSR readGroup x position reverse table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x position reverse table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_position_reverse = new TBQSR_cellPositionRev*[numReadGroups];

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	maxPos = 0;
	int p;
	std::string tmpF;
	std::getline(file, tmpF); //skip header

	//parse file to get max position
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//get quality
			p = stringToInt(vec[1]);
			if(p > maxPos) maxPos = p;
		}
	}

	//create corresponding objects and object to check if we will initialize all positions!
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		BQSR_cells_readGroup_position_reverse[i] = new TBQSR_cellPositionRev[maxPos];
		isListed[i] = new bool[maxPos];
		for(int p=0; p<maxPos; ++p){
			if(considerPosition) BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
			else BQSR_cells_readGroup_position_reverse[i][p].init(BQSR_cells_readGroup_quality, qualityIndex);
			isListed[i][p] = false;
		}
	}

	//rewind file to beginning
	file.clear();
	file.seekg(0, file.beg); //rewind file to beginning
	std::getline(file, tmpF); //skip header
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				p = stringToInt(vec[1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_position_reverse[readGroup][p-1].set(alpha, vec[4]);
				isListed[readGroup][p-1] = true;
			}
		}
	}

	//check if we miss positions
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int p=0; p<maxPos; ++p){
			if(!isListed[i][p]) throw "Position " + toString(p+1) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRPositionReverse")) positionReverseConverged = false;
	else positionReverseConverged = true;
	considerPositionReverse = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering positions reverse up to " + toString(maxPos));
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTable(TParameters & params){
	if(params.parameterExists("BQSRContext")) initializeBQSRReadGroupContextTableFromFile(params);
	else {
		contextConverged = false;
		considerContext = params.parameterExists("estimateBQSRContext");
		if(considerContext){
			estimatetionRequired = true;
			logfile->list("Considering context");
			BQSR_cells_readGroup_context = new TBQSR_cellContext*[numReadGroups];
			for(int r=0; r<numReadGroups; ++r){
				BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
				for(int c=0; c<numContexts; ++c){
					if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex);
					else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
					else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex);
					else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex);
				}
			}
		} else {
			BQSR_cells_readGroup_context = NULL;
		}
	}
}

void TRecalibrationBQSR::initializeBQSRReadGroupContextTableFromFile(TParameters & params){
	std::string filename = params.getParameterString("BQSRContext");
	logfile->listFlush("Constructing BQSR readGroup x context table from file '" + filename + "' ...");
	std::ifstream file(filename.c_str());
	if(!file) throw "Failed to open BQSR readGroup x context table from file '" + filename + "'!";

	//construct for each read group in bam file
	BQSR_cells_readGroup_context = new TBQSR_cellContext*[numReadGroups];
	for(int r=0; r<numReadGroups; ++r){
		BQSR_cells_readGroup_context[r] = new TBQSR_cellContext[numContexts];
		for(int c=0; c<numContexts; ++c){
			if(considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, BQSR_cells_readGroup_position_reverse, qualityIndex);
			else if(considerPosition && !considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position, qualityIndex);
			else if(!considerPosition && considerPositionReverse) BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, BQSR_cells_readGroup_position_reverse, qualityIndex);
			else BQSR_cells_readGroup_context[r][c].init(BQSR_cells_readGroup_quality, qualityIndex);
		}
	}

	//create object to check of all contexts have been initialized!
	bool** isListed = new bool*[numReadGroups];
	for(int i=0; i<numReadGroups; ++i){
		isListed[i] = new bool[numContexts];
		for(int c=0; c<numContexts; ++c){
			isListed[i][c] = false;
		}
	}

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	std::string tmpF;
	std::getline(file, tmpF); //skip header
	int context;
	double alpha;
	int readGroup;

	//now parse file again and set empirical quality
	while(file.good() && !file.eof()){
		++lineNum;
		fillVectorFromLineWhiteSpaceSkipEmpty(file, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 5) throw "Found " + toString(vec.size()) + " instead of 5 columns in '" + filename + "' on line " + toString(lineNum) + "!";
			//set quality and empirical error rate
			readGroup = findReadGroupIndex(vec[0]);
			if(readGroup >= 0){ //returns -1 if read group does not exist
				context = genoMap.getContext(vec[1][0], vec[1][1]);
				alpha = stringToDouble(vec[3]);
				BQSR_cells_readGroup_context[readGroup][context].set(alpha, vec[4]);
				isListed[readGroup][context] = true;
			}
		}
	}

	//check if we miss contexts
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int c=0; c<numContexts; ++c){
			if(!isListed[i][c]) throw "Context " + genoMap.getContextString(c) + " is not listed for read group '" + it->ID + "' in file '" + filename + "'!";
		}
		delete[] isListed[i];
	}
	delete[] isListed;

	//set that no estimation is not required, unless requested
	if(params.parameterExists("estimateBQSRContext")) contextConverged = false;
	else contextConverged = true;
	considerContext = true;

	//done!
	logfile->write(" done!");
	logfile->conclude("Considering context");
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
		}
		else if(considerPositionReverse && !positionReverseConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				if((*it)->posInReadRev >= maxPos) throw "Position of base is > maxPos specified!";
				BQSR_cells_readGroup_position_reverse[(*it)->readGroup][(*it)->posInReadRev].addBase(*it, refBase);
			}
		} else if(considerContext && !contextConverged){
			for(std::vector<TBase*>::iterator it = site.bases.begin(); it != site.bases.end(); ++it){
				BQSR_cells_readGroup_context[(*it)->readGroup][(*it)->context].addBase(*it, refBase);
			}
		}
	}
}

bool TRecalibrationBQSR::estimateEpsilon(){
	//estimate epsilon for quality and readgroup, if not yet done
	estimationConverged = true;
	int numCellsNotConverged = 0;
	double maxF = 0.0;

	//readGroup x quality
	//-------------------------------------------------------
	if(!qualityConverged){
		logfile->listFlush("Estimating epsilon for readGroup x quality table ...");
		for(int i=0; i<numReadGroups; ++i){
			for(int j=0; j<qualityIndex->numQ; ++j){
				if(!BQSR_cells_readGroup_quality[i][j].estimate(convergenceThreshold, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_quality[i][j].F > maxF) maxF = BQSR_cells_readGroup_quality[i][j].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			qualityConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			qualityConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * qualityIndex->numQ));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!qualityConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int j=0; j<qualityIndex->numQ; ++j){
					BQSR_cells_readGroup_quality[i][j].empty();
				}
			}
			estimationConverged = false;
		} else {
			if(!considerPosition && !considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position, if not yet done
	//-------------------------------------------------------
	if(considerPosition && !positionConverged){
		logfile->listFlush("Estimating epsilon for readGroup x position table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position[i][p].estimate(convergenceThreshold, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_position[i][p].F > maxF) maxF = BQSR_cells_readGroup_position[i][p].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			positionConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position[i][p].empty();
				}
			}
			estimationConverged = false;
		} else {
			if(!considerPositionReverse && !considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for position reverse, if not yet done
	//-------------------------------------------------------
	if(considerPositionReverse && !positionReverseConverged){
		logfile->listFlush("Estimating epsilon for readGroup x position reverse table ...");
		numCellsNotConverged = 0;

		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				if(!BQSR_cells_readGroup_position_reverse[i][p].estimate(convergenceThreshold, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_position_reverse[i][p].F > maxF) maxF = BQSR_cells_readGroup_position_reverse[i][p].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			positionReverseConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			positionReverseConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * maxPos));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!positionReverseConverged){
			//empty all cells
			for(int i=0; i<numReadGroups; ++i){
				for(int p=0; p<maxPos; ++p){
					BQSR_cells_readGroup_position_reverse[i][p].empty();
				}
			}
			estimationConverged = false;
		} else {
			if(!considerContext) estimationConverged = true;
			else estimationConverged = false;
		}
		return estimationConverged;
	}

	//estimate epsilon for context
	//-------------------------------------------------------
	if(considerContext && !contextConverged){
		logfile->listFlush("Estimating epsilon for quality x context table ...");
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				if(!BQSR_cells_readGroup_context[r][c].estimate(convergenceThreshold, minObservations)){
					++numCellsNotConverged;
					if(BQSR_cells_readGroup_context[r][c].F > maxF) maxF = BQSR_cells_readGroup_context[r][c].F;
				}
			}
		}

		//report
		logfile->write(" done!");
		if(numCellsNotConverged == 0){
			contextConverged = true;
			logfile->list("Estimation converged in all cells!");
		} else {
			contextConverged = false;
			int percent = 100.0 * ((double) numCellsNotConverged / (double) (numReadGroups * numContexts));
			logfile->conclude("Estimation has not yet converged in " + toString(numCellsNotConverged) + " cells (" + toString(percent) + "%)");
		}
		logfile->conclude("Largest F = " + toString(maxF));

		//set status
		if(!contextConverged){
			//empty all cells
			for(int r=0; r<numReadGroups; ++r){
				for(int c=0; c<numContexts; ++c){
					BQSR_cells_readGroup_context[r][c].empty();
				}
			}
			estimationConverged = false;
		} else {
			estimationConverged = true;
		}
		return estimationConverged;
	}

	//return true on final convergence
	return estimationConverged;
}

void TRecalibrationBQSR::writeToFile(std::string filenameTag){
	//write readGroup x Quality table
	std::string filename = filenameTag + "_BQSR_ReadGroup_Quality_Table.txt";
	logfile->listFlush("Writing BQSR readGroup x quality table to '" + filename + "' ...");
	std::ofstream out(filename.c_str());
	if(!out) throw "Failed to open file '" + filename + "' for writing!";
	out << "ReadGroup\tQualityScore\tEventType\tEmpiricalQuality\tObservations";
	out << "\tfirstDerivative\tSecondDerivative\tF\thasConverged";
	out << "\n";
	BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
	for(int i=0; i<numReadGroups; ++i, ++it){
		for(int q=0; q<qualityIndex->numQ; ++q){
			out << it->ID << "\t" << qualityIndex->getQuality(q) << "\tM\t" << makePhred(BQSR_cells_readGroup_quality[i][q].curEstimate) << "\t" << BQSR_cells_readGroup_quality[i][q].getNumObsForPrinting();
			//for debugging: also print derivatives, F and whether is has converged
			out << "\t" << BQSR_cells_readGroup_quality[i][q].firstDerivativeSave << "\t" << BQSR_cells_readGroup_quality[i][q].secondDerivativeSave << "\t" << BQSR_cells_readGroup_quality[i][q].F << "\t" << BQSR_cells_readGroup_quality[i][q].estimationConverged;
			out << "\n";
		}
	}
	out.close();
	logfile->write(" done!");

	//write readGroup x position table
	if(considerPosition){
		filename = filenameTag + "_BQSR_ReadGroup_Position_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x position table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tPosition\tEventType\tScaling\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int i=0; i<numReadGroups; ++i, ++it){
			for(int p=0; p<maxPos; ++p){
				out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position[i][p].curEstimate << "\t" << BQSR_cells_readGroup_position[i][p].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_position[i][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position[i][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position[i][p].F << "\t" << BQSR_cells_readGroup_position[i][p].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}

	//write readGroup x position_rev table
	if(considerPositionReverse){
		filename = filenameTag + "_BQSR_ReadGroup_Position_Reverse_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x position reverse table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tPosition\tEventType\tScaling\tObservations\n";
		BamTools::SamReadGroupIterator it = bamHeader->ReadGroups.Begin();
		for(int i=0; i<numReadGroups; ++i, ++it){
			for(int p=0; p<maxPos; ++p){
				out << it->ID << "\t" << p+1 << "\tM\t" << BQSR_cells_readGroup_position_reverse[i][p].curEstimate << "\t" << BQSR_cells_readGroup_position_reverse[i][p].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_position_reverse[i][p].firstDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[i][p].secondDerivativeSave << "\t" << BQSR_cells_readGroup_position_reverse[i][p].F << "\t" << BQSR_cells_readGroup_position_reverse[i][p].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}

	//write readGroup x context table
	if(considerContext){
		filename = filenameTag + "_BQSR_ReadGroup_Context_Table.txt";
		logfile->listFlush("Writing BQSR readGroup x text table to '" + filename + "' ...");
		std::ofstream out(filename.c_str());
		if(!out) throw "Failed to open file '" + filename + "' for writing!";
		out << "ReadGroup\tContext\tEventType\tScaling\tObservations\n";
		it = bamHeader->ReadGroups.Begin();
		for(int r=0; r<numReadGroups; ++r, ++it){
			for(int c=0; c<numContexts; ++c){
				out << it->ID << "\t" << genoMap.getContextString(c) << "\tM\t" << BQSR_cells_readGroup_context[r][c].curEstimate << "\t" << BQSR_cells_readGroup_context[r][c].getNumObsForPrinting();
				//for debugging: also print derivatives, F and whether is has converged
				out << "\t" << BQSR_cells_readGroup_context[r][c].firstDerivativeSave << "\t" << BQSR_cells_readGroup_context[r][c].secondDerivativeSave << "\t" << BQSR_cells_readGroup_context[r][c].F << "\t" << BQSR_cells_readGroup_context[r][c].estimationConverged;
				out << "\n";
			}
		}
		out.close();
		logfile->write(" done!");
	}
}

bool TRecalibrationBQSR::allConverged(){
	if(!qualityConverged) return false;
	if(considerPosition && !positionConverged) return false;
	if(considerPositionReverse && !positionReverseConverged) return false;
	if(considerContext && !contextConverged) return false;
	return true;
}

void TRecalibrationBQSR::reopenEstimation(){
	//resets all cells not to have converged
	for(int i=0; i<numReadGroups; ++i){
		for(int q=0; q<qualityIndex->numQ; ++q){
			BQSR_cells_readGroup_quality[i][q].reopenEstimation();
		}
	}
	qualityConverged = false;

	//also for position
	if(considerPosition){
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position[i][p].reopenEstimation();
			}
		}
	}
	positionConverged = false;

	//reverse position
	if(considerPositionReverse){
		for(int i=0; i<numReadGroups; ++i){
			for(int p=0; p<maxPos; ++p){
				BQSR_cells_readGroup_position_reverse[i][p].reopenEstimation();
			}
		}
	}
	positionReverseConverged = false;

	//and context
	if(considerContext){
		for(int r=0; r<numReadGroups; ++r){
			for(int c=0; c<numContexts; ++c){
				BQSR_cells_readGroup_context[r][c].reopenEstimation();
			}
		}
	}
	contextConverged = false;
}

double TRecalibrationBQSR::getErrorRate(TBase* base){
	double q = BQSR_cells_readGroup_quality[base->readGroup][qualityIndex->getIndex(base->quality)].curEstimate;
	if(considerPosition) q *= BQSR_cells_readGroup_position[base->readGroup][base->posInRead].curEstimate;
	if(considerPositionReverse) q *= BQSR_cells_readGroup_position_reverse[base->readGroup][base->posInReadRev].curEstimate;
	if(considerContext) q *= BQSR_cells_readGroup_context[base->readGroup][base->context].curEstimate;
	if(q > 1.0) q = 1.0; //make sure the scaling does not lead to errors > 1.0!
	return q;
}

int TRecalibrationBQSR::getQuality(TBase* base){
	double q = getErrorRate(base);
	//transform to quality
	return makePhredInt(q);
}

