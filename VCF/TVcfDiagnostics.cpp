/*
 * TTVcfDiagnostics.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */

#include "TVcfDiagnostics.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"


namespace VCF {

using coretools::str::toString;
using coretools::instances::parameters;
using coretools::instances::logfile;

//---------------------------------------------------
// TAnnotator
//---------------------------------------------------

TVcfDiagnostics::TVcfDiagnostics() {
	chr = -1;

	// open vcf file
	std::string vcfFileName = parameters().get("vcf");
	auto ending      = coretools::str::readAfterLast(vcfFileName, '.');
	bool isZipped           = (ending == "gz");
	_openVCF(vcfFileName, isZipped);

	// read output name
	auto tmp = coretools::str::readBeforeLast(vcfFileName, ".vcf");
	_outName        = parameters().get("out", tmp);
	logfile().list("Writing output files with prefix '" + _outName + "'. (specify with 'out')");
}

void TVcfDiagnostics::_openVCF(const std::string &VCFName, bool isZipped) {
	// open vcf file
	if (isZipped) {
		logfile().list("Reading vcf from gzipped file '" + VCFName + "'.");
	} else {
		logfile().list("Reading vcf from file '" + VCFName + "'.");
	}

	_vcfFile.openStream(VCFName, isZipped);

	// enable parsers
	_vcfFile.enablePositionParsing();
	_vcfFile.enableVariantParsing();
	_vcfFile.enableVariantQualityParsing();
	_vcfFile.enableFormatParsing();
	_vcfFile.enableSampleParsing();
	_vcfFile.enableInfoParsing();
}

void TVcfDiagnostics::vcfToInvariantBed() {
	// TODO: check what this function is supposed to do; switch to TBed!
	// open vcf file
	logfile().list("Writing sites that are invariant across individuals to bed:");

	// open output
	coretools::TOutputFile bedFile(_outName + ".bed.gz");
	logfile().list("Writing files to '", bedFile.name(), "'");

	// parse vcf file
	std::array<char, 4> bases   = {'A', 'C', 'G', 'T'};
	int counter                 = 0;
	int curStartRegion          = -1;
	bool previousStartIsVariant = false;
	std::string curChr;
	bool updateStart = true;
	long lastPosition;
	while (_vcfFile.next()) {
		++counter;
		if (updateStart) {
			curStartRegion = _vcfFile.position();
			curChr         = _vcfFile.chr();
			updateStart    = false;
		}
		if (std::find(bases.begin(), bases.end(), _vcfFile.getRefAllele()[0]) == bases.end()) {
			continue; // ignore indels
		}

		if (curChr != _vcfFile.chr()) {
			// is previous site invariant? -> write to file
			if (!previousStartIsVariant) { bedFile.writeln(curChr, curStartRegion - 1, lastPosition); }
			// update start directly, don't just set to true
			curStartRegion = _vcfFile.position();
			curChr         = _vcfFile.chr();
		}

		genometools::Base allele = _vcfFile.getFirstAlleleOfSample(0);
		for (int i = 0; i < _vcfFile.numSamples(); ++i) {
			if (_vcfFile.getFirstAlleleOfSample(i) != allele || _vcfFile.getSecondAlleleOfSample(i) != allele) {
				// there was a variant site, is previous site invariant? -> write to file
				if (!previousStartIsVariant && counter != 1) {
					bedFile.writeln(_vcfFile.chr(), curStartRegion - 1, _vcfFile.positionZeroBased());
				}
				updateStart            = true;
				previousStartIsVariant = true;
				break;
			} else {
				previousStartIsVariant = false;
			}
		}

		// in case of chr change, this is needed to write last region of previous chr
		lastPosition = _vcfFile.position();
	}

	// write last region to file
	if (!previousStartIsVariant) {
		bedFile.writeln(_vcfFile.chr(), curStartRegion - 1, _vcfFile.positionZeroBased());
	}
	bedFile.close();
}

int TVcfDiagnostics::findLastPassedFilterIndex(int obsValue, std::vector<int> &filtersAscendingOrder) {
	int lastPassedFilterIndex = 0;
	for (unsigned int i = 0; i < filtersAscendingOrder.size(); ++i) {
		if (obsValue < filtersAscendingOrder[i])
			break;
		else
			lastPassedFilterIndex = i;
	}
	return lastPassedFilterIndex;
}

void TVcfDiagnostics::assessAllelicImbalance() {
	// output
	if (!_vcfFile.formatColExists("AD")) {
		logfile().list("VCF File ", _vcfFile.filename, " does not have allelic depth field!");
		return;
	}
	logfile().list("Writing files to '" + _outName + "_allelicDepth.txt'");

	// limit input?
	int maxDP = parameters().get<int>("maxDepth", 100);
	logfile().list("Ignoring sites with depth larger than " + toString(maxDP) + ".");

	int inputLines = parameters().get<int>("inputLines", -1);
	if (inputLines <= 0) {
		logfile().list("Reading whole vcf.");
	} else
		logfile().list("Limiting input to " + toString(inputLines) + " lines.");

	// initialize tables
	logfile().startIndent("Initializing count tables:");
	std::string qualityString = parameters().get<std::string>("qualities", "0,10,20,30,40,50");
	std::vector<int> qualities;
	coretools::str::fillContainerFromString(qualityString, qualities, ',');

	std::vector<TCountTable *> countTables;
	for (unsigned int i = 0; i < qualities.size(); ++i) {
		// countTables.emplace_back(
		// new TCountTable(maxDP, maxDP, _outName + "_qual" + toString(qualities[i]) + "_allelicDepth.txt", logfile));
	}
	logfile().endIndent();

	// temp variables
	int counter = 0;
	int numHet  = 0;

	logfile().startIndent("Parsing vcf file:");
	while (_vcfFile.next()) {
		++counter;

		// get allelic depth at hetero sites across all samples. One table per quality filter.
		for (int i = 0; i < _vcfFile.numSamples(); ++i) {
			if (!_vcfFile.sampleIsMissing(i) && _vcfFile.sampleIsHeteroRefNonref(i)) {
				++numHet;

				// which position in table does site correspond to?
				std::vector<std::string> tmp;
				std::string tag = "AD";
				coretools::str::fillContainerFromString(_vcfFile.getSampleContentAt(tag, i), tmp, ',');
				int numRef = coretools::str::fromString<int>(tmp[0]);
				int numAlt = coretools::str::fromString<int>(tmp[1]);
				if (numRef == 0 || numAlt == 0)
					UERROR("Call at position ", _vcfFile.position(),
						   " is heterozygous but reference or alternative allelic depth is 0!");
				if (_vcfFile.depthAsIntNoCheckForMissingSample("DP", i) > maxDP) {
					logfile().warning("DP is " + toString(_vcfFile.depthAsIntNoCheckForMissingSample("DP", i)) +
					                  " at pos " + toString(_vcfFile.position()) + ". This site will be ignored.");
					continue;
				}

				// add count to correct table
				int quality = coretools::str::fromString<int>(_vcfFile.getSampleContentAt("GQ", i));
				int index   = findLastPassedFilterIndex(quality, qualities);
				for (int i = 0; i < (index + 1); ++i) { ++(countTables.at(i))->table[numRef][numAlt]; }
			}
		}
		if (counter % 1000000 == 0)
			logfile().list("read " + toString(counter) + " lines, " + toString(numHet) +
			               " heterozygous sites were found.");
		if (inputLines != -1 && counter > inputLines) {
			logfile().list("reached input limit!");
			break;
		}
	}

	logfile().endIndent("done!");

	std::string description = "ref_depth/alt_depth";
	std::string rowPrefix   = "ref_";
	std::string colPrefix   = "alt_";

	// write tables
	for (unsigned int i = 0; i < countTables.size(); ++i) {
		countTables[i]->writeTable(description, rowPrefix, colPrefix);
		delete countTables[i];
	}

	// clean up
	//	for(int i=0; i<countTables.size(); ++i){
	//		delete countTables.at(i);
	//	}
}

void TVcfDiagnostics::fixIntAsFloat() {
	// open vcf file
	logfile().list("Fixing integers that are printed as floats:");


	// open output file
	std::string filename = _outName + (std::string) "_fixed.vcf.gz";

	if (!_vcfFile.formatColExists("GP")) {
		logfile().list("VCF File ", _vcfFile.filename, " does not have a GP field!");
		logfile().list("Will just copy the file");
		if (!std::filesystem::copy_file(_vcfFile.filename, filename))
			UERROR("Failed to copy '", _vcfFile.filename, "' to '", filename, "'!");
		return;
	}

	gz::ogzstream out(filename.c_str());
	if (!out) UERROR("Failed to open outputfile '", filename, "'!");
	_vcfFile.setOutStream(out);
	_vcfFile.writeHeaderVCF_4_0();

	// tmp vars
	int counter = 0;
	std::vector<int> vec;

	// parse VCF

	logfile().startIndent("Parsing vcf file:");
	while (_vcfFile.next()) {
		++counter;

		// fix GP field
		std::string gp = _vcfFile.fieldContentAsString("GP", 0);
		coretools::str::fillContainerFromString(gp, vec, ',');
		gp = std::to_string(vec[0]) + ',' + std::to_string(vec[1]) + ',' + std::to_string(vec[2]);
		_vcfFile.updateField("GP", gp, 0);

		// write output
		_vcfFile.writeLine();

		// report progress
		if (counter % 1000000 == 0) logfile().list("read " + toString(counter) + " lines.");
	}
	logfile().endIndent();
};

}; // namespace VCF
