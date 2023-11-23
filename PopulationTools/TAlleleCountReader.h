/*
 * TAlleleCountReader.h
 *
 *  Created on: Jan 24, 2023
 *      Author: raphael
 */

#ifndef TALLELECOUNTREADER_H_
#define TALLELECOUNTREADER_H_

#include <string>
#include <vector>

#include "coretools/Files/TFile.h"

#include "genometools/GenotypeTypes.h"

namespace PopulationTools{

struct alleleCounts{ 
    uint32_t minor;
    uint32_t total;
};

struct alleleCountVector{
    std::string chr;
    uint64_t pos;
    genometools::Base majorAllele{genometools::Base::N};
    genometools::Base minorAllele{genometools::Base::N};    
    std::vector<alleleCounts> _alleleCounts;

    alleleCounts& operator[](size_t index){ return _alleleCounts[index]; }
    const alleleCounts& operator[](size_t index) const { return _alleleCounts[index]; }

    size_t numPopulations() const { return _alleleCounts.size(); }
};

//-------------------------------------------------
//
//-------------------------------------------------
class TAlleleCountReader{
private:
    coretools::TInputFile _file;
        
    std::vector<std::string> _populationNames;
    alleleCountVector _alleleCountVec;
    bool _hasAlleles = false;
    int _firstPopulationColumn;
    
public:
    using value_type      = alleleCountVector;
	using const_reference = const alleleCountVector&;

    TAlleleCountReader() = default;
    TAlleleCountReader(std::string_view filename);
    void open(std::string_view filename);       
    void close();

    // getters    
    bool empty(){ return !_file.isOpen(); }
    bool hasAlleles(){ return _hasAlleles; }
    size_t numPopulations(){ return _populationNames.size(); }
    const std::vector<std::string>& populationNames(){
        return _populationNames;
    };
    const_reference front(){
        return _alleleCountVec;    
    }

    //read next
    void popFront();

    //iterators
    using iterator = coretools::TLazyIterator<TAlleleCountReader>;
    iterator begin(){ return iterator{this}; }
    iterator end(){ return iterator{}; }
};



}; //end namespace

#endif /* TALLELECOUNTREADER_H_ */
