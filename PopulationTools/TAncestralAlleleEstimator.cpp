/*
 * TGenotypeFreqencies.cpp
 *
 *  Created on: Jan 26, 2023
 *      Author: raphael
 */

#include "TAncestralAlleleEstimator.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TFastaWriter.h"
#include "coretools/Main/TLog.h"

namespace PopulationTools{

using coretools::instances::parameters;
using coretools::instances::logfile;

//--------------------------
// TAncestralAlleleEstimator
//--------------------------
TAncestralAlleleEstimator::TAncestralAlleleEstimator(){
    
    _minorCountMax = parameters().get<size_t>("minorCountMaximum", 0);
    logfile().list("Setting maximum count of minor allele to ", _minorCountMax, ". (parameter 'minorCountMaximum')");
    _totalCountMin = parameters().get<size_t>("totalCountMinimum", 0);
    logfile().list("Setting minimum total allele count to ", _totalCountMin, ". (parameter 'totalCountMinimum')");

    const auto alleleCountsFileName = parameters().get<std::string>("alleleCounts");

    logfile().list("Opening alleleCounts file ", alleleCountsFileName, ". (parameter 'alleleCounts')");
    _alleleCounts.open(alleleCountsFileName);

    filename = coretools::str::readBeforeLast(alleleCountsFileName, "_alleleCounts.txt.gz");

    const auto fastaIndexFileName = parameters().get<std::string>("fastaIndex", filename + ".fasta.fai");
    logfile().list("Opening FASTA index file ", fastaIndexFileName, ". (parameter 'fastaIndex')");
    _fastaIndex.parse(fastaIndexFileName);
}

void TAncestralAlleleEstimator::run(){

    const auto outputFileName = parameters().get<std::string>("out", filename + "_ancestralAlleles.fasta");
    logfile().list("Writing ancestral alleles to file '", outputFileName, "'. (parameter 'out')");

    genometools::TFastaWriter writer(outputFileName, _fastaIndex.lineBases(0));

    std::vector<std::string> populationNames = _alleleCounts.populationNames();
    const auto population = parameters().get<std::string>("population", populationNames[0]);
    //calculating index of population in vector
    auto it = std::find(populationNames.begin(), populationNames.end(), population);
    size_t index;
    if (it != populationNames.end()){
        index = it - populationNames.begin();
    } else {
        throw coretools::TUserError("Population '", population, "' was not found in alleleCounts file!");
    }
    logfile().list("Using population '", populationNames[index],"' for estimation of ancestral allele. (parameter 'population')");

    for (auto &chr: _fastaIndex.chromosomes()){
        writer.newContig(chr.name());
        for (size_t pos = 1; pos <= chr.length(); pos++){
            if(pos == _alleleCounts.front().pos && chr.name() == _alleleCounts.front().chr){
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
