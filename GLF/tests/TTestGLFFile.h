//
// Created by caduffm on 9/1/20.
//

#ifndef ATLAS_TTESTGLFFILE_H
#define ATLAS_TTESTGLFFILE_H

#include <stdint.h>
#include <string>
#include <vector>

#include "TChromosomes.h"
#include "TGLF.h"
#include "TGenomePosition.h"
#include "TGenotypeData.h"

namespace GLF {

//--------------------------------------
// TTestGLFFile
//--------------------------------------
class TTestGLFFile{
protected:
    //header info
    std::string _header;

    //GLF file for writing
    std::string _filename;
    GLF::TGlfWriter _glfFile;

    // chromosomes
    genometools::TChromosomes _chromosomes;
    std::vector<genometools::TChromosome>::iterator _dummyCurChr;

    //tmp vars for dummy sites
    long _dummyPos; long _dummyDist; long _dummyMaxDist;
    uint32_t _dummyDepth; uint32_t _dummyMaxDepth;
    uint8_t _dummy_RMS_mappingQual; uint32_t _dummyMax_RMS_mappingQual;
    GenotypeLikelihoods::TGenotypeLikelihoods _dummyGenotypeLikelihoods;
    GenotypeLikelihoods::TGenotypeLikelihoodsHaploid _dummyGenotypeLikelihoodsHaploid;

    void _iteratePosition();
    void _iterateDepth();
    void _iterateRMSMappingQual();
    virtual void _iterateGenotypeLikelihoods(uint32_t curDepth);

    //storage of written alignments
    std::vector<genometools::TGenomePosition> _writtenPositions;
    std::vector<uint32_t> _writtenDepths;
    std::vector<uint8_t> _writtenRMSMappingQualities;
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _writtenGenotypeLikelihoods;
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _writtenGenotypeLikelihoodsWithMissingSites;

    void _initialize(const std::vector<uint32_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy);

public:
    TTestGLFFile(){};

    void openOutput(const std::string & Filename, std::vector<uint32_t>& ChrLength);
    void openOutput(const std::string & Filename, std::vector<uint32_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy);
    void closeOutput();

    void writeSite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods, uint8_t RMS_mappingQuality);
    void writeNewChromosome();
    // write dummy sites with shuffling
    void writeDummySites(uint32_t numSites);
    void writeDummySite(long pos);
    void writeDummySite(long pos, uint32_t depth);
    void writeDummySite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

    //getters
    std::string filename()const { return _filename; };
    genometools::TChromosomes& chromosomes(){ return _chromosomes; };
    std::vector<genometools::TChromosome>::iterator curChr(){return _dummyCurChr; };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator beginGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.begin(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator endGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.end(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator beginGenotypeLikelihoodsWithMissingSites(){ return _writtenGenotypeLikelihoodsWithMissingSites.begin(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator endGenotypeLikelihoodsWithMissingSites(){ return _writtenGenotypeLikelihoodsWithMissingSites.end(); };
    GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoodsWithMissingSites(uint32_t index){ return _writtenGenotypeLikelihoodsWithMissingSites[index]; };
    std::vector<genometools::TGenomePosition>::iterator beginPositions(){ return _writtenPositions.begin(); };
    std::vector<genometools::TGenomePosition>::iterator endPositions(){ return _writtenPositions.end(); };
    genometools::TGenomePosition position(uint32_t index){ return _writtenPositions[index]; };
    std::vector<uint32_t>::iterator beginDepths(){ return _writtenDepths.begin(); };
    std::vector<uint32_t>::iterator endDepths(){ return _writtenDepths.end(); };
    uint32_t depth(uint32_t index){ return _writtenDepths[index]; };
    std::vector<uint8_t>::iterator beginRMSMappingQualities(){ return _writtenRMSMappingQualities.begin(); };
    std::vector<uint8_t>::iterator endRMSMappingQualities(){ return _writtenRMSMappingQualities.end(); };
    uint8_t RMSMappingQuality(uint32_t index){ return _writtenRMSMappingQualities[index]; };
};

}; // end namespace

#endif //ATLAS_TTESTGLFFILE_H
