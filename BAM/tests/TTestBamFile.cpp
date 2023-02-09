/*
 * TTestBamFile.cpp
 *
 *  Created on: Aug 10, 2020
 *      Author: wegmannd
 */

#include "TTestBamFile.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>

#include "TCigar.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "coretools/Main/globalConstants.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Types/strongTypes.h"

namespace BAM{

//--------------------------------------
// TTestBamFile
//--------------------------------------
TTestBamFile::TTestBamFile(const std::vector<uint32_t> ChrLength, uint32_t NumReadGroups){
	_initialize(ChrLength, NumReadGroups);
};

TTestBamFile::TTestBamFile(const std::string & Filename, const std::vector<uint32_t> ChrLength, uint32_t NumReadGroups){
	_initialize(ChrLength, NumReadGroups);

	//open BAM file for writing
	_filename = Filename;
	_bamFile.open(_filename, _header, _chromosomes, _readGroups);
};

void TTestBamFile::_initialize(const std::vector<uint32_t> ChrLength, uint32_t NumReadGroups){
	using genometools::Base;
	_header.set("1.6", "coordinate", "none", "none");
	_initializeChromosomes(ChrLength);
	_initializeReadGroups(NumReadGroups);

	// tmp vars
	_dummySequence = {Base::A, Base::A, Base::A, Base::C, Base::C, Base::C, Base::G, Base::G, Base::G, Base::T, Base::T,
			  Base::T, Base::A, Base::C, Base::G, Base::T, Base::T, Base::G, Base::C, Base::A, Base::A, Base::A,
			  Base::C, Base::G, Base::T, Base::G, Base::G, Base::C, Base::C, Base::G, Base::T, Base::G, Base::A,
			  Base::C, Base::A, Base::C, Base::C, Base::G, Base::T, Base::C, Base::G, Base::A, Base::C, Base::A,
			  Base::G, Base::G, Base::T, Base::G, Base::C, Base::C, Base::A, Base::C, Base::A, Base::C, Base::A,
			  Base::G, Base::T, Base::G, Base::G, Base::C, Base::A, Base::A, Base::A, Base::T, Base::T, Base::G,
			  Base::G, Base::C, Base::C, Base::G, Base::G, Base::T, Base::G, Base::C, Base::A, Base::A, Base::A,
			  Base::C, Base::C, Base::A, Base::A, Base::A, Base::C, Base::C, Base::A, Base::A, Base::G, Base::G,
	                  Base::T, Base::T, Base::G, Base::C, Base::C, Base::C, Base::G};
	_dummySequenceStart = _dummySequence.begin();
	for(auto p = genometools::BaseQuality::min(); p < genometools::BaseQuality::max(); ++p){
		_dummyQualities.emplace_back(genometools::PhredIntProbability(p));
	}
	_dummyQualitiesStart = _dummyQualities.begin();
    _dummyMapQual = 0;
	_dummyMinLength = 10;
	_dummyMaxLength = 50;
	_dummyLength = _dummyMinLength;
	_dummyReadGroup = 0;
	_dummyIsReverseStrand = false;
	_dummyCigarChars = "MID";
	_dummyCigarPos = 0;
	_dummyFlag = 0;
    _counter = 0;
};

void TTestBamFile::_initializeChromosomes(const std::vector<uint32_t> ChrLength){
	uint32_t refID = 0;
	for(auto& len : ChrLength){
		++refID;
		_chromosomes.appendChromosome("Chr" + coretools::str::toString(refID), len);
	}
};

void TTestBamFile::_initializeReadGroups(uint32_t NumReadGroups){
	for(uint32_t i=0; i<NumReadGroups; ++i){
		BAM::TReadGroup& rg = _readGroups.add("ReadGroup" + coretools::str::toString(i));
		rg.sequencingCenter_CN = coretools::__GLOBAL_APPLICATION_NAME__ + " " + coretools::__GLOBAL_APPLICATION_VERSION__;
		rg.description_DS = "Simulated with commit " + coretools::__GLOBAL_APPLICATION_COMMIT__;
		rg.sample_SM = "TestSample";
		rg.sequencingTechnology_PL = "ILLUMINA";
	}
};

void TTestBamFile::_iterateMappingQuality(){
    _dummyMapQual = (_dummyMapQual + 4) % 255;
}

void TTestBamFile::_iterateReadGroupAndReverseStrand(){
    _dummyReadGroup = (_dummyReadGroup + 1) % _readGroups.size();
    _dummyIsReverseStrand = !_dummyIsReverseStrand;
}

void TTestBamFile::_iterateCigar(BAM::TCigar & cigar, uint32_t length){
    std::string s = _dummyCigarChars.substr(_dummyCigarPos, 3);
    while(s.length() < 3){
        s += _dummyCigarChars.substr(0, 3 - s.length());
    }
    // construct cigar: choose number of S, M, I and D
    cigar.add('S', std::floor(length / 8));
    for (char const &c: s) {
        cigar.add(c, std::floor(length / 3));
    }
    if (length > cigar.lengthRead())
        cigar.add('S',  length - cigar.lengthRead()); // fill rest with S (right)

    // update
    _dummyCigarPos = (_dummyCigarPos + 1) % _dummyCigarChars.size();
}

void TTestBamFile::_iterateLength(){
    _dummyLength += 3;
    if(_dummyLength > _dummyMaxLength){
        _dummyLength = _dummyMinLength + _dummyLength % _dummyMaxLength;
    }
}

void TTestBamFile::_iterateFlags() {
    // sam flags are bits -> can be summed to one value; increasing this value means a different combination of bits
    // there are 12 flags -> 2^12 = 4096 - 1 = 4095 is the last valid combination
    if (_dummyFlag.asInt() < 4095){
        _dummyFlag = _dummyFlag.asInt() + 1;
        while (!_dummyFlag.isValid() && _dummyFlag.isPaired()){ // some sums are invalid -> only accept valid sums. Also, we only want to simulate single-end reads in here
            _dummyFlag = _dummyFlag.asInt() + 1;
            if (_dummyFlag.asInt() >= 4095)
                _dummyFlag = 0;
        }
    } else _dummyFlag = 0;
}

void TTestBamFile::openOutput(const std::string & Filename){
	//open BAM file for writing
	_filename = Filename;
	_bamFile.open(_filename, _header, _chromosomes, _readGroups);
};

void TTestBamFile::closeOutput(){
	_bamFile.close();
};

void TTestBamFile::_storeAlignment(const BAM::TAlignment & alignment){
	//store for later comparisons
	_writtenAlignments.emplace_back(alignment);
};

void TTestBamFile::writeAlignment(const BAM::TAlignment & alignment){
    //write to BAM
    _bamFile.writeAlignment(alignment);
};

BAM::TAlignment TTestBamFile::_constructAlignment(const std::vector<genometools::Base> & sequence, const std::vector<genometools::PhredIntProbability> & qualities, const genometools::TGenomePosition & position, const BAM::TCigar & cigar, uint32_t readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag){
	BAM::TAlignment alignment(position);
    alignment.setName("alignment_" + coretools::str::toString(_counter));
    alignment.setSequenceQualities(cigar, sequence, qualities);
    alignment.setMappingQuality(genometools::PhredIntProbability(_dummyMapQual));
    alignment.setReadGroup(readGroup);
    if (complicatedSamFlag) {
        alignment.setSamFlags(BAM::TSamFlags(_dummyFlag));
        _iterateFlags();
    } else
        alignment.setIsReverseStrand(isReverseStrand);

    return alignment;
}

void TTestBamFile::writeDummyAlignment(const genometools::TGenomePosition & position, const BAM::TCigar & cigar, uint32_t readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag){
	//extract sequence / qualities
	auto s = _constructFrom(_dummySequence, _dummySequenceStart, cigar.lengthRead());
	auto q = _constructFrom(_dummyQualities, _dummyQualitiesStart, cigar.lengthRead());

	// iterate mapping quality
	_iterateMappingQuality();

	// fill alignment
    BAM::TAlignment alignment = _constructAlignment(s, q, position, cigar, readGroup, isReverseStrand, complicatedSamFlag);

    // store and write
    _storeAlignment(alignment);
    writeAlignment(alignment);
    _counter++;
};

void TTestBamFile::writeDummyAlignment(const genometools::TGenomePosition & position, const BAM::TCigar & cigar, const bool & complicatedSamFlag){
	writeDummyAlignment(position, cigar, _dummyReadGroup, _dummyIsReverseStrand, complicatedSamFlag);

    //iterate
    _iterateReadGroupAndReverseStrand();
};

void TTestBamFile::writeDummyAlignment(const genometools::TGenomePosition & position, uint32_t length, const bool & complicatedSamFlag){
    // iterate over order of M, I and D
    BAM::TCigar cigar;
    _iterateCigar(cigar, length);

    // now write
	writeDummyAlignment(position, cigar, complicatedSamFlag);
};

void TTestBamFile::writeDummyAlignment(const genometools::TGenomePosition & position, const bool & complicatedSamFlag){
	writeDummyAlignment(position, _dummyLength, complicatedSamFlag);
	_iterateLength();
};

uint32_t TTestBamFile::_computeDistanceBetweenAlignments(uint32_t numAlignments){
    uint32_t usableLength = _chromosomes.referenceLength() - _chromosomes.size() * _dummyMaxLength;
    uint32_t dist = usableLength / ((double) numAlignments + 1);
    if(dist > _chromosomes.minLength() - _dummyMaxLength){
        dist = _chromosomes.minLength() - _dummyMaxLength;
    }
    return dist;
}

void TTestBamFile::writeDummyAlignments(uint32_t numAlignments, const bool & complicatedSamFlag){
	//get distance between alignments
    uint32_t dist = _computeDistanceBetweenAlignments(numAlignments);

    genometools::TGenomePosition position;
	auto chr = _chromosomes.begin();

	for(uint32_t i=0; i<numAlignments; ++i){
		//iterate position
		position += dist;
		if(position + _dummyLength > chr->chrEnd){
			++chr;
			if(chr == _chromosomes.end()){
				std::cout << "ERROR A" << std::endl;
				DEVERROR("chromosome reached end!");
			}

			position = chr->chrStart;
		}
		writeDummyAlignment(position, complicatedSamFlag);
	}
}

void TTestBamFile::writeDummyAlignment(const genometools::Base &oneBase, const genometools::PhredIntProbability &oneQual, const genometools::TGenomePosition &position,
                                  const BAM::TCigar &cigar, uint32_t readGroup,
                                  const bool &isReverseStrand) {
    // expand oneBase to string
	std::vector<genometools::Base> seq(cigar.lengthRead(), oneBase);
	std::vector<genometools::PhredIntProbability> qual(cigar.lengthRead(), oneQual);

    //write alignment
    BAM::TAlignment alignment(position);
    alignment.setSequenceQualities(cigar, seq, qual);
    alignment.setReadGroup(readGroup);
    alignment.setIsReverseStrand(isReverseStrand);
    writeAlignment(alignment);


    /*
    std::cout << "----------------------------------" << std::endl;
    std::cout << "refID = " << position.refID() << std::endl;
    std::cout << "position = " << position.position() << std::endl;
    std::cout << "length = " << cigar.lengthRead() << std::endl;
    std::cout << "seq = " << seq << std::endl;
    std::cout << "qual = " << qual << std::endl;
    std::string cigg;
    for (auto & it : cigar){
        std::string cigger(it.length, it.type);
        cigg += cigger;
    }
    std::cout << "cigar = " << cigg << std::endl;
    std::cout << "readGroup = " << (int) readGroup << std::endl;
    std::cout << "isReverseStrand = " << isReverseStrand << std::endl;
    */
};


void TTestBamFile::writeDummyAlignment(const genometools::Base &oneBase, const genometools::PhredIntProbability &oneQual, const genometools::TGenomePosition & position, const BAM::TCigar & cigar){
    writeDummyAlignment(oneBase, oneQual, position, cigar, _dummyReadGroup, _dummyIsReverseStrand);

    //iterate
    _iterateReadGroupAndReverseStrand();
}

void TTestBamFile::writeDummyAlignment(const genometools::Base &oneBase, const genometools::PhredIntProbability &oneQual, const genometools::TGenomePosition &position,
                                  uint32_t length) {
    // construct cigar
    BAM::TCigar cigar;
    _iterateCigar(cigar, length);

    writeDummyAlignment(oneBase, oneQual, position, cigar);
};

void TTestBamFile::writeDummyAlignment(const genometools::Base &oneBase, const genometools::PhredIntProbability &oneQual, const genometools::TGenomePosition &position) {
    writeDummyAlignment(oneBase, oneQual, position, _dummyLength);

    // iterate length
    _iterateLength();
};

//--------------------------------------
// TTestBamFilePairedEnd
//--------------------------------------

TTestBamFilePairedEnd::TTestBamFilePairedEnd(const std::vector<uint32_t> ChrLength, uint32_t NumReadGroups) : TTestBamFile(ChrLength, NumReadGroups){
    _dummyFlag = 1;
}

TTestBamFilePairedEnd::TTestBamFilePairedEnd(const std::string &Filename, const std::vector<uint32_t> ChrLength,
                                             uint32_t NumReadGroups) : TTestBamFile(Filename, ChrLength, NumReadGroups) {
    _dummyFlag = 1;

}

void TTestBamFilePairedEnd::_iterateFlags() {
    // sam flags are bits -> can be summed to one value; increasing this value means a different combination of bits
    // there are 12 flags -> 2^12 = 4096 - 1 = 4095 is the last valid combination
    if (_dummyFlag.asInt() < 4095){
        _dummyFlag = _dummyFlag.asInt() + 1;
        while (!_dummyFlag.isValid() || !_dummyFlag.isPaired()){ // some sums are invalid -> only accept valid sums. Also, we only want to simulate paired-end reads in here
            _dummyFlag = _dummyFlag.asInt() + 1;
            if (_dummyFlag.asInt() >= 4095)
                _dummyFlag = 0;
        }
    } else _dummyFlag = 1;
}

void TTestBamFilePairedEnd::writeDummyAlignment(const genometools::TGenomePosition & position, const BAM::TCigar & cigar, uint32_t readGroup, const bool & isReverseStrand, const bool & complicatedSamFlag){
    //extract sequence / qualities
	auto s = _constructFrom(_dummySequence, _dummySequenceStart, cigar.lengthRead());
	auto q = _constructFrom(_dummyQualities, _dummyQualitiesStart, cigar.lengthRead());

    // iterate mapping quality
    _iterateMappingQuality();

    BAM::TAlignment alignment = _constructAlignment(s, q, position, cigar, readGroup, isReverseStrand, complicatedSamFlag);

    // store, but DON'T write! We don't know information about mate yet
    _storeAlignment(alignment);
    _counter++;
};

BAM::TAlignment & TTestBamFilePairedEnd::_pickFirstMate(std::vector<bool> & used){
    // pick first mate: first alignment that has not been used yet
    uint32_t indexMate1 = 0;
    for (uint32_t j = 0; j < _writtenAlignments.size(); j++){
        if (!used[j]){
            used[j] = true;
            indexMate1 = j;
            break;
        }
    }
    return _writtenAlignments[indexMate1];
}

BAM::TAlignment & TTestBamFilePairedEnd::_pickSecondMate(uint32_t refIDMate1, std::vector<bool> & used){
    // pick second mate: take last position that is still on same chromosome as mate1 and is not used -> that way, fragment lengths will vary a lot
    uint32_t indexMate2 = 0;
    for (uint32_t j = _writtenAlignments.size() - 1; j >= 0; j--){
        if (_writtenAlignments[j].refID() == refIDMate1 && !used[j]) {
            used[j] = true;
            indexMate2 = j;
            break;
        }
    }
    return _writtenAlignments[indexMate2];
}


void TTestBamFilePairedEnd::writeDummyAlignments(uint32_t numAlignments, const bool & complicatedSamFlag){
    if (numAlignments % 2 != 0)
        DEVERROR("For simplicity, numAlignments should be an even number!");

    // create alignments and store
    BAM::TTestBamFile::writeDummyAlignments(numAlignments, complicatedSamFlag);

    // now match forward and reverse read
    std::vector<bool> used(_writtenAlignments.size(), false);
    for (uint32_t i = 0; i < numAlignments/2.; i++){
        BAM::TAlignment & mate1 = _pickFirstMate(used);
        BAM::TAlignment & mate2 = _pickSecondMate(mate1.refID(), used);

        // set mate information
        mate1.setMateGenomicPosition(mate2);
        mate1.setInsertSize(mate2.position() - mate1.position());
        mate1.setIsReverseStrand(false);
        mate1.setIsRead1(true); mate1.setIsRead2(false);

        mate2.setMateGenomicPosition(mate1);
        mate2.setInsertSize(mate2.position() - mate1.position());
        mate2.setIsReverseStrand(true);
        mate2.setIsRead1(false); mate2.setIsRead2(true);
    }

    // finally: write
    for (auto & alignment : _writtenAlignments){
        writeAlignment(alignment);
    }
}

}; //end namespace
