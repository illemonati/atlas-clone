#include "TAlleleCountReader.h"

namespace PopulationTools{

TAlleleCountReader::TAlleleCountReader(std::string_view filename){
    open(filename); 
}

void TAlleleCountReader::open(std::string_view filename){   
    if(_file.isOpen()){
        DEVERROR("Allele count file is already open!");
    }
    //open file
    _file.open(filename, coretools::FileType::Header);

    // parse header	
    const auto header = _file.header();

    if(header[0] != "chr" || header[1] != "pos"){
        UERROR("Allele count file '", filename, "' lacks columns 'chr' and 'pos' at beginning! Are you providing the correct file?");
    }

    if(header[2] == "ref" && header[3] == "alt"){
        _hasAlleles = true;
        _firstPopulationColumn = 4;
    } else {
        _hasAlleles = false;
        _firstPopulationColumn = 2;
    }

    _populationNames.assign(header.begin() + _firstPopulationColumn, header.end());
    _alleleCountVec._alleleCounts.resize(numPopulations());

    popFront();
}

void TAlleleCountReader::close(){
    _file.close();
}

// read next line, if there is no next line -> close file
void TAlleleCountReader::popFront(){

    //set chr and position
    _alleleCountVec.chr = _file.get(0);
    _file.fill(1, _alleleCountVec.pos);    
    if(_hasAlleles){
        _alleleCountVec.majorAllele = genometools::char2base(_file.get(2)[0]);
        _alleleCountVec.minorAllele = genometools::char2base(_file.get(3)[0]);
    }

    //extract counts
    for(size_t p = 0; p < _populationNames.size(); ++p){    
        coretools::str::TSplitter sp (_file.get(p + _firstPopulationColumn), '/');
        coretools::str::fromString(sp.front(), _alleleCountVec[p].minor);
        sp.popFront();
        coretools::str::fromString(sp.front(), _alleleCountVec[p].total);        
    }    
	_file.popFront();
}

}; //end namespace
