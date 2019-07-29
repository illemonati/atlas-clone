#include "TAtlasTestList.h"

TAtlasTestList::TAtlasTestList(){
    //fill map of tests
    //NOTE: order will be the order in which test will be run if initialized from TParameters
    testMap.emplace_back("empty", &createInstance<TAtlasTest>);
    testMap.emplace_back("pileup", &createInstance<TAtlasTest_pileup>);
    testMap.emplace_back("allelicDepth", &createInstance<TAtlasTest_allelicDepth>);
    testMap.emplace_back("recalSimulation", &createInstance<TAtlasTest_recalSimulation>);
    testMap.emplace_back("BQSRSimulation", &createInstance<TAtlasTest_BQSRSimulation>);
    testMap.emplace_back("qualityTransformationRecalPlain", &createInstance<TAtlasTest_qualityTransformationRecalPlain>);
    testMap.emplace_back("qualityTransformationRecalBinned", &createInstance<TAtlasTest_qualityTransformationRecalBinned>);
    testMap.emplace_back("PMDEmpiric", &createInstance<TAtlasTest_PMDEmpiric>);
    testMap.emplace_back("theta", &createInstance<TAtlasTest_theta>);
    testMap.emplace_back("invariantBed", &createInstance<TAtlasTest_invariantBed>);
    testMap.emplace_back("mergePairs", &createInstance<TAtlasTest_mergePairs>);
    testMap.emplace_back("filter", &createInstance<TAtlasTest_filter>);
    testMap.emplace_back("removeSoftClips", &createInstance<TAtlasTest_removeSoftClips>);
    testMap.emplace_back("librarieCompareFunctions", &createInstance<TAtlasTestLibrarieCompareFunctions>);
    testMap.emplace_back("librarieCompareSpeed", &createInstance<TAtlasTestLibrarieSpeed>);

    //fill map of test suites
    //Note: suites and tests within suites will be initialized in this order!
    testSuites["exampleSuite"] = {"empty", "pileup"};

    //automatically create a test suite "all"
    testSuites["all"] = {};
    for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt)
        testSuites["all"].push_back(testMapIt->first);

}

TAtlasTestList::~TAtlasTestList(){
    for(std::vector<TAtlasTest*>::iterator it=initializedTests.begin(); it!=initializedTests.end(); ++it)
        delete (*it);
}

bool TAtlasTestList::initializeTest(const std::string &name, TParameters &params, TLog *logfile){
    //check if test exists
    if(testInitialized(name))
        return true;

    //if not create it
    for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
        if(testMapIt->first == name){
            initializedTests.push_back( testMapIt->second(params, logfile) );
            return true;
        }
    }

    //only reached if test does not exist
    return false;
}

void TAtlasTestList::parseTests(TParameters &params, TLog *logfile){
    for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
        if(params.parameterExists(testMapIt->first))
            initializeTest(testMapIt->first, params, logfile);
    }
}

bool TAtlasTestList::parseSuites(TParameters &params, TLog *logfile){
    bool returnVal = false;
    for(std::map<std::string, std::vector<std::string> >::iterator it = testSuites.begin(); it != testSuites.end(); ++it){
        if(params.parameterExists(it->first)){
            for(size_t i=0; i<it->second.size(); ++i)
                returnVal &= initializeTest(it->second[i], params, logfile);
        }
    }
    return returnVal;
}

size_t TAtlasTestList::size(){
    return initializedTests.size();
}

bool TAtlasTestList::testExists(const std::string name){
    for(testMapIt = testMap.begin(); testMapIt != testMap.end(); ++testMapIt){
        if(testMapIt->first == name)
            return true;
    }
    return false;
}

bool TAtlasTestList::testInitialized(const std::string name){
    for(testIt=initializedTests.begin(); testIt!=initializedTests.end(); ++testIt){
        if((*testIt)->name() == name)
            return true;
    }
    return false;
}

void TAtlasTestList::printTestToLogfile(TLog *logfile){
    for(testIt=initializedTests.begin(); testIt!=initializedTests.end(); ++testIt){
        logfile->list((*testIt)->name());
    }
}

TAtlasTest *TAtlasTestList::operator[](size_t num){
    if(num < initializedTests.size())
        return initializedTests[num];
    else throw "Test number out of range!";
}
