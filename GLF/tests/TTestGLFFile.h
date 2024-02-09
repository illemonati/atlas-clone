//
// Created by caduffm on 9/1/20.
//

#ifndef ATLAS_TTESTGLFFILE_H
#define ATLAS_TTESTGLFFILE_H

#include <stdint.h>
#include <string>
#include <vector>

#include "genometools/GenomePositions/TChromosomes.h"
#include "TGLF.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "GenotypeData.h"

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
    size_t _dummyDepth; size_t _dummyMaxDepth;
    uint8_t _dummy_RMS_mappingQual; size_t _dummyMax_RMS_mappingQual;
    GenotypeLikelihoods::TGenotypeLikelihoods _dummyGenotypeLikelihoods{1.};
    GenotypeLikelihoods::TGenotypeLikelihoods _dummyGenotypeLikelihoodsHaploid{1.};

    void _iteratePosition();
    void _iterateDepth();
    void _iterateRMSMappingQual();
    virtual void _iterateGenotypeLikelihoods(size_t curDepth);

    //storage of written alignments
    std::vector<genometools::TGenomePosition> _writtenPositions;
    std::vector<size_t> _writtenDepths;
    std::vector<uint8_t> _writtenRMSMappingQualities;
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _writtenGenotypeLikelihoods;
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods> _writtenGenotypeLikelihoodsWithMissingSites;

    void _initialize(const std::vector<size_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy);

public:
    TTestGLFFile(){};

    void openOutput(const std::string & Filename, std::vector<size_t>& ChrLength);
    void openOutput(const std::string & Filename, std::vector<size_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy);
    void closeOutput();

    void writeSite(long pos, size_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods, uint8_t RMS_mappingQuality);
    void writeNewChromosome();
    // write dummy sites with shuffling
    void writeDummySites(size_t numSites);
    void writeDummySite(long pos);
    void writeDummySite(long pos, size_t depth);
    void writeDummySite(long pos, size_t depth, GenotypeLikelihoods::TGenotypeLikelihoods & genotypeLikelihoods);

    //getters
    std::string filename()const { return _filename; };
    genometools::TChromosomes& chromosomes(){ return _chromosomes; };
    std::vector<genometools::TChromosome>::iterator curChr(){return _dummyCurChr; };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator beginGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.begin(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator endGenotypeLikelihoods(){ return _writtenGenotypeLikelihoods.end(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator beginGenotypeLikelihoodsWithMissingSites(){ return _writtenGenotypeLikelihoodsWithMissingSites.begin(); };
    std::vector<GenotypeLikelihoods::TGenotypeLikelihoods>::iterator endGenotypeLikelihoodsWithMissingSites(){ return _writtenGenotypeLikelihoodsWithMissingSites.end(); };
    GenotypeLikelihoods::TGenotypeLikelihoods genotypeLikelihoodsWithMissingSites(size_t index){ return _writtenGenotypeLikelihoodsWithMissingSites[index]; };
    std::vector<genometools::TGenomePosition>::iterator beginPositions(){ return _writtenPositions.begin(); };
    std::vector<genometools::TGenomePosition>::iterator endPositions(){ return _writtenPositions.end(); };
    genometools::TGenomePosition position(size_t index){ return _writtenPositions[index]; };
    std::vector<size_t>::iterator beginDepths(){ return _writtenDepths.begin(); };
    std::vector<size_t>::iterator endDepths(){ return _writtenDepths.end(); };
    size_t depth(size_t index){ return _writtenDepths[index]; };
    std::vector<uint8_t>::iterator beginRMSMappingQualities(){ return _writtenRMSMappingQualities.begin(); };
    std::vector<uint8_t>::iterator endRMSMappingQualities(){ return _writtenRMSMappingQualities.end(); };
    uint8_t RMSMappingQuality(size_t index){ return _writtenRMSMappingQualities[index]; };
};

}; // end namespace

#endif //ATLAS_TTESTGLFFILE_H
