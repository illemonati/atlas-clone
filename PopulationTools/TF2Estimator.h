//
// Created by reynac on 2/22/22.
//

#ifndef ATLAS_TF2ESTIMATOR_H
#define ATLAS_TF2ESTIMATOR_H

#include "TParameters.h"
#include "TPopulationLikelihoods.h"
#include "TRandomGenerator.h"
#include "mathFunctions.h"

namespace PopulationTools {
//------------------------------------------------
//TF2Estimator
//------------------------------------------------
    class TF2Estimator {
    private:
        coretools::TLog *_logfile;
        coretools::TRandomGenerator *_randonGenerator;
        std::string _outname;

        //vcf-file
        std::string _vcfFilename;
        VCF::TVcfFileSingleLine _vcfFile;
        bool _limitLines;
        uint64_t _maxNumLines;

        //samples
        TPopulationSamples _samples;

        //genotype data
        //THWPopulations _populations;
    public:
        TF2Estimator(coretools::TParameters &Parameters, coretools::TLog *Logfile, coretools::TRandomGenerator *RandonGenerator);
        void calculateF2();
    };

    //--------------------------------------
// Tasks
//--------------------------------------
    class TTask_calculateF2:public coretools::TTask{
    public:
        TTask_calculateF2(){ _explanation = "Calculate F2 between different samples"; };

        void run(coretools::TParameters & Parameters, coretools::TLog* Logfile){
            TF2Estimator runF2(Parameters, Logfile, _randomGenerator);
            runF2.calculateF2();
        };
    };

}

#endif //ATLAS_TF2ESTIMATOR_H
