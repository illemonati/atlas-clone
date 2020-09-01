//
// Created by caduffm on 9/1/20.
//

#include "TTestGLFFile.h"

namespace TestUtilities {

TTestGLFFile::TTestGLFFile(const std::vector<uint32_t>& ChrLength) {
    _initialize(ChrLength);
}


TTestGLFFile::TTestGLFFile(const std::string &Filename, const std::vector<uint32_t>& ChrLength) {
    _initialize(ChrLength);

    openOutput(Filename);
}

void TTestGLFFile::_initialize(const std::vector<uint32_t>& ChrLength) {
    uint32_t refID = 0;
    for(auto& len : ChrLength){
        ++refID;
        _chromosomes.appendChromosome("Chr" + toString(refID), len);
    }
    _dummyCurChr = _chromosomes.begin();

    _dummyDepth = 0; _dummyMaxDepth = 50;
    _dummyPos = 0; _dummyDist = 0; _dummyMaxDist = 0;
    _dummy_RMS_mappingQual = 0; _dummyMax_RMS_mappingQual = 20;
}

void TTestGLFFile::openOutput(const std::string & Filename){
    //open GLF file for writing
    _filename = Filename;
    _header = "";
    _glfFile.open(_filename, _header);
}

void TTestGLFFile::closeOutput() {
    _glfFile.close();
}

void TTestGLFFile::_iterateDepth() {
    _dummyDepth = (_dummyDepth + 1) % _dummyMaxDepth;
}

void TTestGLFFile::_iterateRMSMappingQual() {
    _dummy_RMS_mappingQual = (_dummy_RMS_mappingQual + 3) % _dummyMax_RMS_mappingQual;
}

void TTestGLFFile::_iteratePosition() {
    _dummyDist = (_dummyDist + 7) % _dummyMaxDist;
    _dummyPos += _dummyDist + 1; // + 1 because % can return 0, but then we would like distance of 1

    // next chromosome?
    if (BAM::TGenomePosition(_dummyCurChr->refID(), _dummyPos) > _dummyCurChr->chrEnd){
        ++_dummyCurChr;
        if(_dummyCurChr == _chromosomes.end()){
            throw std::runtime_error("void TTestBamFile::writeDummyAlignments(const uint32_t & numAlignments): chromosome reached end!");
        }

        _dummyPos = _dummyCurChr->chrStart.position();
        writeNewChromosome(*_dummyCurChr);
    }
}

void TTestGLFFile::_iterateGenotypeLikelihoods(uint32_t curDepth) {
    std::vector<Base> possibleBases = {A, C, G, T};
    uint8_t indexPossibleBases = 0;

    double error = 0.01;
    double maxError = 0.25;

    // simulate a site with different base and errors
    std::vector<GenotypeLikelihoods::TBaseData> bases;
    for (uint32_t b = 0; b < curDepth; b++){
        // pick a true base
        Base trueBase = possibleBases[indexPossibleBases];
        // simulate a base with error
        bases.emplace_back(trueBase, error);

        // iterate, such that next base gets different values
        indexPossibleBases = (indexPossibleBases + 3) % possibleBases.size();
        error += 0.02;
        if (error > maxError)
            error = 0.01;
    }

    // fill genotype likelihood
    _dummyGenotypeLikelihoods.fill(bases);
}

void TTestGLFFile::writeDummySites(const uint32_t &numSites) {
    // compute maximal distance between sites (such that positions never exceed the last position of the last chromosome)
    uint32_t usableLength = _chromosomes.referenceLength();
    if (numSites > usableLength)
        throw std::runtime_error("Too many sizes (" + toString(numSites) + ") for chromosomes of total length " + toString(usableLength) + "!");
    _dummyMaxDist = usableLength / numSites;
    if(_dummyMaxDist > _chromosomes.minLength() - 1){
        _dummyMaxDist = _chromosomes.minLength() - 1;
    }

    // write first chromosome
    writeNewChromosome(*_dummyCurChr);

    // go over all sites, simulate some data and write
    for(uint32_t i=0; i<numSites; ++i){
        _iteratePosition();
        writeDummySite(_dummyPos);
    }
}

void TTestGLFFile::writeDummySite(long pos) {
    writeDummySite(pos, _dummyDepth);

    _iterateDepth();
}

void TTestGLFFile::writeDummySite(long pos, uint32_t depth) {
    _iterateGenotypeLikelihoods(depth);

    writeDummySite(pos, depth, _dummyGenotypeLikelihoods);
}

void TTestGLFFile::writeDummySite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods){
    writeSite(pos, depth, genotypeLikelihoods, _dummy_RMS_mappingQual);

    _iterateRMSMappingQual();
}

void TTestGLFFile::writeSite(long pos, uint32_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods, uint8_t RMS_mappingQual) {
    // write to glf
    _glfFile.writeSite(pos, depth, RMS_mappingQual, genotypeLikelihoods);

    // ... and store, for later comparisons
    _writtenPositions.push_back(pos);
    _writtenDepths.push_back(depth);
    _writtenGenotypeLikelihoods.push_back(genotypeLikelihoods);
    _writtenRMSMappingQualities.push_back(RMS_mappingQual);
}

void TTestGLFFile::writeNewChromosome(const BAM::TChromosome & chromosome) {
    _glfFile.newChromosome(chromosome);
}














}; // end namespace