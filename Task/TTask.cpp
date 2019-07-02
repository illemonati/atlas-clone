#include "TTask.h"

void TTask::initializeRandomGenerator(TParameters &parameters, TLog *logfile){
    if(!randomGeneratorInitialized){
        logfile->listFlush("Initializing random generator ...");

        if(parameters.parameterExists("fixedSeed")){
            randomGenerator=new TRandomGenerator(parameters.getParameterLong("fixedSeed"), true);
        } else if(parameters.parameterExists("addToSeed")){
            randomGenerator=new TRandomGenerator(parameters.getParameterLong("addToSeed"), false);
        } else {
            randomGenerator=new TRandomGenerator();
        }
        logfile->write(" done with seed " + toString(randomGenerator->usedSeed) + "!");
        randomGeneratorInitialized = true;
    }
}

TTask::TTask(){
    randomGenerator = nullptr;
    randomGeneratorInitialized = false;
}

TTask::~TTask(){
    if(randomGeneratorInitialized){
        delete randomGenerator;
    }
}

void TTask::run(std::string taskName, TParameters &parameters, TLog *logfile){
    logfile->startIndent(_explanation + " (task = " + taskName + "):");

    //print citations
    printCitations(logfile);

    initializeRandomGenerator(parameters, logfile);

    //now run task
    run(parameters, logfile);
    logfile->endIndent();
}

void TTask::run(TParameters &parameters, TLog *logfile){
    throw "Base task can not be run!";
}

std::string TTask::explanation(){ return _explanation; }

void TTask::printCitations(TLog *logfile){
    //write citations, if there are any
    if(!_citations.empty()){
        logfile->startIndent("Relevant citations:");
        for(std::vector<std::string>::iterator it = _citations.begin(); it != _citations.end(); ++it){
            logfile->list(*it);
        }
        logfile->endIndent();
    }
}
