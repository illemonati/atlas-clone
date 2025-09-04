/*
 * TCompareVCF.cpp
 *
 *  Created on: May 7, 2019
 *      Author: wegmannd
 */

#include "TVcfCompare.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace VCF{

using genometools::Base;
using genometools::Genotype;
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::str::toString;
using coretools::user_assert;

//--------------------------------------------------------------
// TGenotypeComparisonTable
//--------------------------------------------------------------
//add haploid genotypes
void TGenotypeComparisonTable::add(const Base b1, const Base b2){
	++_counts[coretools::index(b1)][coretools::index(b2)];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Base b){
	if(sample == 0){
		++_counts[coretools::index(b)].back();
	} else {
		++_counts.back()[coretools::index(b)];
	}
};

//add diploid genotypes
void TGenotypeComparisonTable::add(Genotype g1, Genotype g2){
	++_counts[_firstDiploidIndex + coretools::index(g1)][_firstDiploidIndex + coretools::index(g2)];
};

void TGenotypeComparisonTable::addOtherMissing(const int sample, const Genotype g){
	if(sample == 0){
		++_counts[_firstDiploidIndex + coretools::index(g)].back();
	} else {
		++_counts.back()[_firstDiploidIndex + coretools::index(g)];
	}
};

//add haploid / diploid combination of genotypes
void TGenotypeComparisonTable::add(const Genotype g1, const Base b2){
	++_counts[_firstDiploidIndex + coretools::index(g1)][coretools::index(b2)];
};

void TGenotypeComparisonTable::add(const Base b1, const Genotype g2){
	++_counts[coretools::index(b1)][_firstDiploidIndex + coretools::index(g2)];
};

//write
void TGenotypeComparisonTable::write(const std::string filename){
	//open output file
	coretools::TOutputFile out(filename);

	//write header
	std::vector<std::string> header = {"vcf1/vcf2"};
	//haploid bases
	for(auto b = Base::min; b < Base::max; ++b){
		header.push_back(toString(b));
	}

	//diploid genotypes
	for(auto g = Genotype::min; g < Genotype::max; ++g){
		header.push_back(toString(g));
	}

	//missing
	header.push_back("N/NN");

	//write header
	out.writeHeader(header);

	//write counts
	for(int i=0; i<_size; i++){
		//write row name
		out << header[i+1];

		for(int j=0; j<_size; j++){
			out << _counts[i][j];
		}

		out.endln();
	};
};

//--------------------------------------------------------------
// TVCFComapreVCF
//--------------------------------------------------------------

TVcfComapreVCF::TVcfComapreVCF(std::string_view filename, std::string_view sampleName){
	//open vcf file
	if(filename.find(".gz") == std::string::npos){
		logfile().list("Reading sample '", sampleName, "' from VCF file '", filename, "'.");
		_vcfFile = std::make_unique<genometools::TVcfFileSingleLine>(filename, false);
	} else {
		logfile().list("Reading sample '", sampleName, "' from gzipped VCF file '", filename, "'.");
		_vcfFile = std::make_unique<genometools::TVcfFileSingleLine>(filename, true);
	}

	_vcfFile->enablePositionParsing();
	_vcfFile->enableFormatParsing();
	_vcfFile->enableSampleParsing();
	_vcfFile->enableVariantParsing();

	//add sample index
	_sampleIndex = _vcfFile->sampleNumber(sampleName);

	//move to first position
	_vcfFile->next();

	coretools::user_assert(!_vcfFile->eof, "Vcf file '", filename, "' is empty!");

	//store first chr
	_parsedChromosomes.push_back(_vcfFile->chr());

	//set filters to zero
	_minDepth = 0;
	_minQual = 0.0;
};


void TVcfComapreVCF::next(){
	_vcfFile->next();
	if(!_vcfFile->eof){
		//store chr if new
		if(_vcfFile->chr() != _parsedChromosomes.back())
			_parsedChromosomes.push_back(_vcfFile->chr());

		//filter
		if(!_vcfFile->sampleIsMissing(_sampleIndex) && _minDepth > 0){
			if(_vcfFile->sampleDepth(_sampleIndex) < _minDepth){
				_vcfFile->setSampleMissing(_sampleIndex);
			}
		}

		if(!_vcfFile->sampleIsMissing(_sampleIndex) && _minQual > 0){
			if(_vcfFile->sampleGenotypeQuality(_sampleIndex) < _minQual){
				_vcfFile->setSampleMissing(_sampleIndex);
			}
		}
	}
};

void TVcfComapreVCF::setFilters(const int MinDepth, const double MinQual){
	_minDepth = MinDepth;
	_minQual = MinQual;
};

bool TVcfComapreVCF::chrParsed(const std::string chr){
	return std::find(_parsedChromosomes.begin(), _parsedChromosomes.end(), chr) != _parsedChromosomes.end();
};

//--------------------------------------------------------------
// TVCFCompare
//--------------------------------------------------------------

void TVcfCompare::addToOtherMissing(TGenotypeComparisonTable & counts, const int sample){
	if(!_vcfFiles[sample].isMissing()){
		if(_vcfFiles[sample].isDiploid()){
			counts.addOtherMissing(sample, _vcfFiles[sample].genotype());
		} else {
			counts.addOtherMissing(sample, _vcfFiles[sample].base());
		}
	}
};

TVcfCompare::TVcfCompare() {
	//open vcf files
	logfile().startIndent("Open VCF files to compare:");
	const auto fileNames   = parameters().get<std::vector<std::string>>("vcf");
	const auto sampleNames = parameters().get<std::vector<std::string>>("samples");

	//currently only implemented for comparing two VCFs
	if (fileNames.size() == 1) {
		user_assert(sampleNames.size() == 2,
							   "VCF comparison requires either one file and two sample names (not ", sampleNames.size(),
							   "), or two files and one sample name!");
	} else if (fileNames.size() == 2) {
		const auto sampleSize = sampleNames.size();
		user_assert(sampleSize == 1 || sampleSize == 2, "VCF comparison requires either one file and two sample names (not ", sampleSize, "), or two files and one sample name!");
	}
	// back is either same as front or the 2nd one, depending whether size = 1 or 2
	_vcfFiles.emplace_back(fileNames.front(), sampleNames.front());
	_vcfFiles.emplace_back(fileNames.back(), sampleNames.back());

	logfile().endIndent();

	//are filters in place?
	int minDepth = parameters().get<int>("minDepth", 0);
	if(minDepth > 0){
		logfile().list("Will consider genotypes with depth < " + toString(minDepth) + " as missing.");
	}
	double minQual = parameters().get("minQual", 0.0);
	if(minQual > 0){
		logfile().list("Will consider genotypes with quality < " + toString(minQual) + " as missing.");
	}

	//limitLines
	if(parameters().exists("limitLines")){
		_limitLines = true;
		logfile().list("Will stop reading after ", _limitLines, " lines.");
		_lineLimit = parameters().get<int>("limitLines");
	}

	//set filters in VCF files
	for(TVcfComapreVCF& it : _vcfFiles){
		it.setFilters(minDepth, minQual);
	}

	_outName = parameters().get("out", "");
	if(_outName.empty()){
		//guess from filename
		//get base name of first VCF file
		_outName = fileNames[0];
		_outName = coretools::str::extractBeforeLast(_outName, '.');
		if(fileNames[0].find(".gz") != std::string::npos){
			//if zipped there is extra .gz
			_outName = coretools::str::extractBeforeLast(_outName, ".");
		}

		if (fileNames.size() > 1) {
			// get base name of second VCF file
			std::string tmp = fileNames[1];
			tmp             = coretools::str::extractBeforeLast(tmp, '.');
			if (fileNames[1].find(".gz") != std::string::npos) {
				// if zipped there is extra .gz
				tmp = coretools::str::extractBeforeLast(tmp, ".");
			}

			_outName += "_" + tmp;
		}
	}

	_outName += "_CallComparison.txt";
}

void TVcfCompare::run(){

	//prepare count table
	TGenotypeComparisonTable counts;

	//now parse VCf files and compare calls
	logfile().startIndent("Parsing vcf file:");
	uint32_t numLines = 0;
	struct timeval start;
	gettimeofday(&start, NULL);
	while((!_vcfFiles[0].eof() || !_vcfFiles[1].eof()) && (!_limitLines  || numLines < _lineLimit)){
		//is one end of file?
		if(_vcfFiles[0].eof()){
			addToOtherMissing(counts, 1);
			_vcfFiles[1].next();
		} else if(_vcfFiles[1].eof()){
			addToOtherMissing(counts, 0);
			_vcfFiles[0].next();
		}

		//are we on the same chromosome?
		else if(_vcfFiles[0].chr() == _vcfFiles[1].chr()){
			//same chromosome
			if(_vcfFiles[0].position() == _vcfFiles[1].position()){
				if(_vcfFiles[0].isMissing()){
					//do not add comparisons where both are missing
					if(!_vcfFiles[1].isMissing()){
						addToOtherMissing(counts, 1);
					}
				} else {
					if(_vcfFiles[1].isMissing()){
						addToOtherMissing(counts, 0);
					} else {
						//both have calls
						if(_vcfFiles[0].isDiploid()){
							if(_vcfFiles[1].isDiploid()){
								counts.add(_vcfFiles[0].genotype(), _vcfFiles[1].genotype());
							} else {
								counts.add(_vcfFiles[0].genotype(), _vcfFiles[1].base());
							}
						} else {
							if(_vcfFiles[1].isDiploid()){
								counts.add(_vcfFiles[0].base(), _vcfFiles[1].genotype());
							} else {
								counts.add(_vcfFiles[0].base(), _vcfFiles[1].base());
							}
						}
					}
				}

				//advance both
				_vcfFiles[0].next();
				_vcfFiles[1].next();
			} else {
				if(_vcfFiles[0].position() < _vcfFiles[1].position()){
					//position is missing in vcfFile1
					addToOtherMissing(counts, 0);
					_vcfFiles[0].next();
				} else {
					//position is missing in vcfFile0
					addToOtherMissing(counts, 1);
					_vcfFiles[1].next();
				}
			}
		} else {
			//we are on different chromosomes
			//has VCFFile1 already parsed chromosome of vcfFile0?
			if(_vcfFiles[1].chrParsed(_vcfFiles[0].chr())){
				//vcfFile1 is on next chromosome, hence, position is missing in vcfFile1
				addToOtherMissing(counts, 0);
				_vcfFiles[0].next();
			} else if(_vcfFiles[0].chrParsed(_vcfFiles[1].chr())){
				//vcfFile0 is on next chromosome, hence, position is missing in vcfFile0
				addToOtherMissing(counts, 1);
				_vcfFiles[1].next();
			} else {
				throw coretools::TUserError("Chromosomes differ between the two VCF files!");
			}
		}

		//report progress
		++numLines;
		if(numLines % 10000 == 0){
			struct timeval end;
			gettimeofday(&end, NULL);
			float runtime = (end.tv_sec  - start.tv_sec)/60.0;
			logfile().list("Parsed " + toString(numLines) + " lines in " + toString(runtime) + " min.");
		}
	}
	struct timeval end;
	gettimeofday(&end, NULL);
	float runtime = (end.tv_sec  - start.tv_sec)/60.0;
	logfile().list("Parsed " + toString(numLines) + " lines in " + toString(runtime) + " min.");
	logfile().list("Reached end of files.");
	logfile().endIndent();

	//write output file

	//writing output file
	logfile().listFlush("Writing counts to file '" + _outName + "' ...");
	counts.write(_outName);
	logfile().done();
};

}; //end namespace
