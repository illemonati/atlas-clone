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
        std::string _outname;

        bool _initialized;
        std::string _vcfFilename;
        VCF::TVcfFileSingleLine _vcfFile;
        bool _limitLines;
        uint64_t _maxNumLines;
        uint32_t _minDepth;
        uint32_t _minNumSamplesWithData;
        uint32_t _minVariantQuality;

        //samples
        TPopulationSamples _samples;
        void _init();
        void _initialize();
        void _openVCF();
        int _filterOnDepth();
        bool _filterSite();
        void _readVCF();

    public:
        TF2Estimator();
        void run();
        void calculateF2(const std::vector<uint64_t>& countsDiff);
        void writeF2(const std::vector<uint64_t>& countsDiff, const std::vector<double>& sampleF2, const std::vector<double>& popF2);
    };

    //--------------------------------------
// Tasks
//--------------------------------------
	class TTask_calculateF2 : public coretools::TTask {
	public:
		TTask_calculateF2() {
			_explanation = "Calculate F2 between different samples, and within and between populations";
		};

		void run() {
			using namespace coretools::instances;
			TF2Estimator F2;
            F2.run();
		};
	};
	} // namespace PopulationTools

#endif //ATLAS_TF2ESTIMATOR_H
