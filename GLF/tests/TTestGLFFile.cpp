//
// Created by caduffm on 9/1/20.
//

#include "TTestGLFFile.h"
#include "GenotypeFunctions.h"

namespace GLF {

void TTestGLFFile::_initialize(const std::vector<size_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy) {
    if (ChrLength.size() != ChrPloidy.size())
        DEVERROR("Size of ChrLength should match size of ChrPloidy!");

    _chromosomes.clear();
    size_t refID = 1;
    for (size_t chr = 0; chr < ChrLength.size(); chr++, refID++){
        _chromosomes.appendChromosome("Chr" + coretools::str::toString(refID), ChrLength[chr], ChrPloidy[chr]);
    }
    _dummyCurChr = _chromosomes.begin() - 1;

    _dummyDepth = 0; _dummyMaxDepth = 50;
    _dummyPos = 0; _dummyDist = 0; _dummyMaxDist = 0;
    _dummy_RMS_mappingQual = 0; _dummyMax_RMS_mappingQual = 20;

    // initialize all entries of _writtenGenotypeLikelihoodsWithMissingSites with missing
    for (size_t pos = 0; pos < _chromosomes.referenceLength(); pos++){
	    _writtenGenotypeLikelihoodsWithMissingSites.emplace_back(1.);
    }
};

void TTestGLFFile::openOutput(const std::string & Filename, std::vector<size_t>& ChrLength){
	std::vector<uint8_t> chrPloidy(ChrLength.size(), 2); // all diploid
	openOutput(Filename, ChrLength, chrPloidy);
};

void TTestGLFFile::openOutput(const std::string & Filename, std::vector<size_t>& ChrLength, const std::vector<uint8_t>& ChrPloidy){
	_initialize(ChrLength, ChrPloidy);

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
    if (genometools::TGenomePosition(_dummyCurChr->refID(), _dummyPos) >= _dummyCurChr->to()){
        writeNewChromosome();
    }
}

void TTestGLFFile::_iterateGenotypeLikelihoods(size_t curDepth) {
	using namespace genometools;
    std::vector<genometools::Base> possibleBases = {Base::A, Base::C, Base::G, Base::T};
    uint8_t indexPossibleBases = 0;

    coretools::Probability error = 0.01;
    coretools::Probability maxError = 0.25;

    // simulate a site with different base and errors
    std::vector<GenotypeLikelihoods::TBaseLikelihoods> bases;
    for (size_t b = 0; b < curDepth; b++){
        // pick a true base
    	genometools::Base trueBase = possibleBases[indexPossibleBases];
        // simulate a base with error
        bases.push_back(GenotypeLikelihoods::fromError(trueBase, error));

        // iterate, such that next base gets different values
        indexPossibleBases = (indexPossibleBases + 3) % possibleBases.size();
        error = error + 0.02;
        if (error > maxError)
            error = 0.01;
    }

    // fill genotype likelihood
    _dummyGenotypeLikelihoods = GenotypeLikelihoods::getGLH(bases);
    /*if (_dummyCurChr->ploidy == 1)
        _dummyGenotypeLikelihoodsHaploid.fill(bases, bases.size());
	else*/
};

void TTestGLFFile::writeDummySites(size_t numSites) {
    // compute maximal distance between sites (such that positions never exceed the last position of the last chromosome)
    size_t usableLength = _chromosomes.referenceLength();
    if (numSites > usableLength)
        DEVERROR("Too many sizes (", numSites, ") for chromosomes of total length ", usableLength, "!");
    _dummyMaxDist = usableLength / numSites;
    if(_dummyMaxDist > _chromosomes.minLength() - 1){
        _dummyMaxDist = _chromosomes.minLength() - 1;
    }

    // write first site
    writeNewChromosome();
    writeDummySite(_dummyPos);

    // go over all sites, simulate some data and write
    for(size_t i=1; i<numSites; ++i){
        _iteratePosition();
        writeDummySite(_dummyPos);
    }
}

void TTestGLFFile::writeDummySite(long pos) {
    writeDummySite(pos, _dummyDepth);

    _iterateDepth();
}

void TTestGLFFile::writeDummySite(long pos, size_t depth) {
    _iterateGenotypeLikelihoods(depth);

    if (_dummyCurChr->ploidy() == 1)
        writeDummySite(pos, depth, _dummyGenotypeLikelihoodsHaploid);
    else
        writeDummySite(pos, depth, _dummyGenotypeLikelihoods);
}

void TTestGLFFile::writeDummySite(long pos, size_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods){
    writeSite(pos, depth, genotypeLikelihoods, _dummy_RMS_mappingQual);

    _iterateRMSMappingQual();
}

void TTestGLFFile::writeSite(long pos, size_t depth, GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods, uint8_t RMS_mappingQual) {
    if (pos < 0 || pos >= _dummyCurChr->length())
        DEVERROR("Position ", pos, " is outside interval [0, lengthChr)!");
    // write to glf
    _glfFile.writeSite(pos, depth, RMS_mappingQual, genotypeLikelihoods);

    long positionInGenome = pos;
    for (auto & chr : _chromosomes){
        if (chr.refID() == _dummyCurChr->refID()) break;
        positionInGenome += chr.length();
    }

    // ... and store, for later comparisons
    _writtenPositions.emplace_back(_dummyCurChr->refID(), pos);
    _writtenDepths.push_back(depth);
    _writtenGenotypeLikelihoods.push_back(genotypeLikelihoods);
    _writtenGenotypeLikelihoodsWithMissingSites[positionInGenome] = genotypeLikelihoods;
    _writtenRMSMappingQualities.push_back(RMS_mappingQual);
}

void TTestGLFFile::writeNewChromosome() {
    ++_dummyCurChr;
    if(_dummyCurChr == _chromosomes.end()){
        DEVERROR("void TTestBamFile::writeDummyAlignments(size_t numAlignments): chromosome reached end!");
    }

    _dummyPos = _dummyCurChr->from().position();

    _glfFile.newChromosome(*_dummyCurChr);
}

}; // end namespace
