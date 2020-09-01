//
// Created by caduffm on 9/1/20.
//

#ifndef ATLAS_TTESTGLFFILE_H
#define ATLAS_TTESTGLFFILE_H

#include "TGLF.h"
#include "globalConstants.h"
#include "stringFunctions.h"

namespace TestUtilities {

//--------------------------------------
// TTestGLFFile
//--------------------------------------
class TTestGLFFile{
protected:
    //header info
    std::string _header;

    // TODO: why is GLF not inside a namespace?

    //GLF file for writing
    std::string _filename;
    TGlfWriter _glfFile;

    // chromosomes
    BAM::TChromosomes _chromosomes;
    std::vector<BAM::TChromosome>::iterator _dummyCurChr;

    //tmp vars for dummy sites
    long _dummyPos; long _dummyDist; long _dummyMaxDist;
    uint32_t _dummyDepth; uint32_t _dummyMaxDepth;
    uint8_t _dummy_RMS_mappingQual; uint32_t _dummyMax_RMS_mappingQual;
    GenotypeLikelihoods::TGenotypeLikelihoods _dummyGenotypeLikelihoods;

    void _iteratePosition();
    void _iterateDepth();
    void _iterateRMSMappingQual();
    virtual void _iterateGenotypeLikelihoods(uint32_t curDepth);

    //storage of written alignments
    std::vector<long> _writtenPositions;
    std::vector<uint32_t> _writtenDepths;
    std::vector<uint8_t> _writtenRMSMappingQualities;
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _writtenGenotypeLikelihoods;

    void _initialize(const std::vector<uint32_t>& ChrLength);

public:
    explicit TTestGLFFile(const std::vector<uint32_t>& ChrLength);
    TTestGLFFile(const std::string & Filename, const std::vector<uint32_t>& ChrLength);

    void openOutput(const std::string & Filename);
    void closeOutput();
    void writeSite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods, uint8_t RMS_mappingQuality);
    void writeNewChromosome(const BAM::TChromosome & chromosome);
    // write dummy sites with shuffling
    void writeDummySites(const uint32_t & numSites);
    void writeDummySite(long pos);
    void writeDummySite(long pos, uint32_t depth);
    void writeDummySite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

    //getters
    std::string filename()const { return _filename; };
    BAM::TChromosomes& chromosomes(){ return _chromosomes; };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator beginGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.begin(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator endGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.end(); };
    std::vector<long>::iterator beginPositions(){ return _writtenPositions.begin(); };
    std::vector<long>::iterator endPositions(){ return _writtenPositions.end(); };
    std::vector<uint32_t>::iterator beginDepths(){ return _writtenDepths.begin(); };
    std::vector<uint32_t>::iterator endDepths(){ return _writtenDepths.end(); };
    std::vector<uint8_t>::iterator beginRMSMappingQualities(){ return _writtenRMSMappingQualities.begin(); };
    std::vector<uint8_t>::iterator endRMSMappingQualities(){ return _writtenRMSMappingQualities.end(); };
};

}; // end namespace

#endif //ATLAS_TTESTGLFFILE_H
