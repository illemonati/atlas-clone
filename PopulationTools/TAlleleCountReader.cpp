#include "TAlleleCountReader.h"

namespace PopulationTools{

TAlleleCountReader::TAlleleCountReader(std::string_view filename){
    open(filename); 
}

void TAlleleCountReader::open(std::string_view filename){
    if(!_reader->isOpen()){
        _chrName.clear();
        _pos = 0;
        _minorAlleleCounts.clear();
        _totalAlleleCounts.clear();
        _populationNames.clear();
        
    }
    _reader = std::make_unique<coretools::TStdReader>(filename);
    // read header
	std::vector<std::string> lineVec;
	file.read(lineVec);
	for (unsigned int i = 2; i < lineVec.size(); ++i) { _populationNames.push_back(lineVec[i]); }
}

// implement function that reads next line, saves different values in its properties, skips empty lines and returns false when eof
bool TAlleleCountReader::readNextLine(){
    std::vector<std::string> lineVec;
	file.read(lineVec);
    _chrName = lineVec[0];
    _pos = coretools::str::convertString<size_t>(lineVec[1]);

    for (unsigned int p = 2; p < lineVec.size(); p++) {
        std::vector<std::string> counts;
        coretools::str::fillContainerFromStringAny(lineVec[p], counts, "/");
        _minorAlleleCounts.push_back(counts[0]);
        _totalAlleleCounts.push_back(counts[1]);
    }
    
}

}; //end namespace