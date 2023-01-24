/*
 * TAlleleCountReader.h
 *
 *  Created on: Jan 24, 2023
 *      Author: raphael
 */

#ifndef TALLELECOUNTREADER_H_
#define TALLELECOUNTREADER_H_

#include <stdint.h>
#include <string>
#include <vector>
#include "coretools/Files/TFile.h"
#include "coretools/Files/TReader.h"

namespace PopulationTools{
    //-------------------------------------------------
class TAlleleCountReader{
private:
    mutable std::unique_ptr<coretools::TReader> _reader = std::make_unique<coretools::TNoReader>();
    std::string _chrName;
    size_t _pos;
    std::vector<std::string> _minorAlleleCounts;
    std::vector<std::string> _totalAlleleCounts;
    std::vector<std::string> _populationNames;
    //put pointer to TInputFile
public:
    TAlleleCountReader() = default;
    TAlleleCountReader(std::string_view filename);
    void open(std::string_view filename);
    bool readNextLine();
    
    void close() { _reader = std::make_unique<coretools::TNoReader>(); }

    // getters
    std::string chrName(){return _chrName;};
    size_t pos(){return _pos;};
    std::vector<std::string> minorAlleleCounts(){return _minorAlleleCounts;};
    std::vector<std::string> totalAlleleCounts(){return _totalAlleleCounts;};
    std::vector<std::string> populationNames(){return _populationNames;};
};
}; //end namespace

#endif /* TALLELECOUNTREADER_H_ */