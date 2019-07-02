#include "TGenotypePrior.h"

//TGenotypePrior
TGenotypePrior::TGenotypePrior(){
    for(int g=0; g<10; ++g)
        genotypePrior[g] = 1.0;
}

TGenotypePrior::~TGenotypePrior(){}

void TGenotypePrior::update(TWindow *window, TLog *logfile){}

double *TGenotypePrior::getPointerToPrior(){ return genotypePrior; }


//TGenotypePriorUniform
TGenotypePriorUniform::TGenotypePriorUniform(){
    for(int g=0; g<10; ++g)
        genotypePrior[g] = 1.0 / 10.0;
}


//TGenotypePriorFixedTheta
TGenotypePriorFixedTheta::TGenotypePriorFixedTheta(double theta, bool EqualBaseFreq, TLog *logfile){
    thetaEstimator = new TThetaEstimator(logfile);
    thetaEstimator->setTheta(theta);
    equalBaseFreq = EqualBaseFreq;
    if(equalBaseFreq){
        TBaseFrequencies freq;
        freq.setEqualBaseFreq();
        thetaEstimator->setBaseFreq(freq);
    }
    thetaEstimator->fillPGenotype(genotypePrior);
}

TGenotypePriorFixedTheta::~TGenotypePriorFixedTheta(){
    delete thetaEstimator;
}

void TGenotypePriorFixedTheta::update(TWindow *window, TLog *logfile){
    if(!equalBaseFreq){
        logfile->listFlush("Estimating base frequencies for prior ...");
        window->estimateBaseFrequencies();
        thetaEstimator->setBaseFreq(window->baseFreq);
        logfile->done();
        thetaEstimator->fillPGenotype(genotypePrior);
    }
}


//TGenotypePriorTheta
void TGenotypePriorTheta::init(TParameters &parameters, std::string &thetaOutputName, TLog *Logfile){
    logfile = Logfile;
    thetaEstimator = new TThetaEstimator(parameters, logfile);
    out.open(thetaOutputName, thetaEstimator, logfile);
}

TGenotypePriorTheta::TGenotypePriorTheta(TParameters &parameters, std::string thetaOutputName, TLog *logfile){
    hasDefaultTheta = false;
    defaultTheta = -1.0;

    init(parameters, thetaOutputName, logfile);
}

TGenotypePriorTheta::TGenotypePriorTheta(TParameters &parameters, std::string thetaOutputName, double DefaultTheta, TLog *logfile){
    hasDefaultTheta = true;
    defaultTheta = DefaultTheta;
    if(defaultTheta < 0.0) throw "Theta must be >= 0.0!";
    init(parameters, thetaOutputName, logfile);	}

TGenotypePriorTheta::~TGenotypePriorTheta(){
    out.close();
    delete thetaEstimator;
}

void TGenotypePriorTheta::update(TWindow *window, TLog *logfile, std::string &chrName){
    logfile->startIndent("Estimating theta and base frequencies:");
    //clear theta estimator
    (*thetaEstimator).clear();

    //adding sites to estimator
    window->addSitesToThetaEstimator(thetaEstimator->pointerToDataContainer());

    //estimate Theta
    if(!thetaEstimator->estimateTheta()){
        if(hasDefaultTheta){
            logfile->conclude("Will use a default theta of " + toString(defaultTheta) + ".");
            thetaEstimator->setTheta(defaultTheta);
            TBaseFrequencies freq;
            freq.setEqualBaseFreq();
            thetaEstimator->setBaseFreq(freq);
        } else
            throw "Please increase window size or provide a default theta!";
    }

    //write results to file
    out.writeWindow(window->chrName, window->start, window->end, thetaEstimator);
    logfile->endIndent();
}
