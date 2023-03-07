/*
 * TGenotypeFreqencies.cpp
 *
 *  Created on: Jan 26, 2023
 *      Author: raphael
 */

#include "TAncestralAlleleEstimator.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TFastaWriter.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Main/TLog.h"

namespace PopulationTools{

using coretools::instances::parameters;
using coretools::instances::logfile;

//--------------------------
// TAncestralAlleleEstimator
//--------------------------
TAncestralAlleleEstimator::TAncestralAlleleEstimator(){
    
    _minorCountMax = parameters().getParameterWithDefault<uint32_t>("minorCountMaximum", 0);
    logfile().list("Setting maximum count of minor allele to ", _minorCountMax, ". (parameter 'minorCountMaximum')");
    _totalCountMin = parameters().getParameterWithDefault<uint32_t>("totalCountMinimum", 0);
    logfile().list("Setting minimum total allele count to ", _totalCountMin, ". (parameter 'totalCountMinimum')");

    const auto alleleCountsFileName = parameters().getParameter<std::string>("alleleCounts");

    logfile().list("Opening alleleCounts file ", alleleCountsFileName, ". (parameter 'alleleCounts')");
    _alleleCounts.open(alleleCountsFileName);

    filename = coretools::str::readBeforeLast(alleleCountsFileName, "_alleleCounts.txt.gz");

    const auto fastaIndexFileName = parameters().getParameterWithDefault<std::string>("fastaIndex", filename + ".fasta.fai");
    logfile().list("Opening FASTA index file ", fastaIndexFileName, ". (parameter 'fastaIndex')");
    _fastaIndex.open(fastaIndexFileName);
}

void TAncestralAlleleEstimator::run(){

    const auto outputFileName = parameters().getParameterWithDefault<std::string>("out", filename + "_ancestralAlleles.fasta");
    logfile().list("Writing ancestral alleles to file '", outputFileName, "'. (parameter 'out')");

    genometools::TFastaWriter writer(outputFileName, _fastaIndex.front().lineBases);

    std::vector<std::string> populationNames = _alleleCounts.populationNames();
    const auto population = parameters().getParameterWithDefault<std::string>("population", populationNames[0]);
    //calculating index of population in vector
    auto it = std::find(populationNames.begin(), populationNames.end(), population);
    uint16_t index;
    if (it != populationNames.end()){
        index = it - populationNames.begin();
    } else {
        UERROR("Population '", population, "' was not found in alleleCounts file!");
    }
    logfile().list("Using population '", populationNames[index],"' for estimation of ancestral allele. (parameter 'population')");

    for (auto &fI: _fastaIndex){
        writer.newContig(fI.contig);
        for (uint64_t pos = 1; pos <= fI.length; pos++){
            if(pos == _alleleCounts.front().pos && fI.contig == _alleleCounts.front().chr){
                if(_alleleCounts.front()[index].minor <= _minorCountMax && _totalCountMin <= _alleleCounts.front()[index].total){
                    writer.write(_alleleCounts.front().majorAllele);
                } else {
                    writer.write(genometools::Base::N);
                }
                _alleleCounts.popFront();
            } else {
                writer.write(genometools::Base::N);
            }
        }
    }
    writer.close();
}

}; //end namespace
