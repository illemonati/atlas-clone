//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "PhredProbabilityTypes.h"
#include "TFile.h"
#include "TLog.h"
#include "TParameters.h"
#include "TPopulation.h"
#include "TPopulationLikelihoods.h"
#include "TSampleLikelihoods.h"
#include "TStorage.h"
#include "TTask.h"

namespace genometools {
class TBed;
}
namespace genometools {
template<typename Type> class TPopulationLikehoodLocus;
}

namespace VCF {

using TSampleLikelihoods = genometools::TSampleLikelihoods<genometools::HighPrecisionPhredIntProbability>;

//------------------------------------------
// TVcfConverter
//------------------------------------------
class TVcfConverter {
protected:
	// filenames
	std::string _vcfName;
	std::string _outName;

	// vcf reader
	genometools::TPopulationLikelihoodReaderLocus _reader;

	// specify population and samples to keep
	genometools::TPopulationSamples _samples;

	void _openVCF();
	void _readSamples();
	void _parseVCF();

	virtual void _initOutputFiles() = 0;
	virtual void _writeHeader();
	virtual void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) = 0;

public:
	TVcfConverter();
	virtual ~TVcfConverter() = default;
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
class TVcfToBeagle : protected TVcfConverter {
private:
	coretools::TOutputFile _beagleFile;
	bool _foundHaploid = false;

	void _writeHeader() override;
	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;
	void _initOutputFiles() override;

	void _writeHaploid(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data, size_t s);
	void _writeDiploid(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data, size_t s);

public:
	TVcfToBeagle();
	~TVcfToBeagle() override = default;

	void run();
};

class TVcfBeagleNew {
public:
	void run();
};

//------------------------------------------
// TVcfToGeno
//------------------------------------------
class TVcfToGeno : protected TVcfConverter {
	// geno format: used in LEA package (https://www.rdocumentation.org/packages/LEA/versions/1.4.0/topics/geno)
	// rows = SNPs, cols = individuals (without delim, just pasted together)
	// Each data point represents the number of copies of reference allele. Missing data encoded by 9.
	// Note: we usually define genotype as number of copies of alternative allele -> might be confusing, I stick
	// to the file format here, i.e. I use number of copies of reference allele.
private:
	coretools::TOutputFile _genoFile;
	coretools::TOutputFile _lociNamesFile;

	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;
	void _initOutputFiles() override;
	void _closeOutputFiles();

public:
	TVcfToGeno();
	~TVcfToGeno() override = default;

	void run();
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------

template<typename T> class TVcfToLFMM : protected TVcfConverter {
protected:
	// genotypes
	coretools::TMultiDimensionalStorage<T, 2> _genotypes;

	// files
	coretools::TOutputFile _lfmmFile;
	coretools::TOutputFile _lociNamesFile;
	std::vector<std::string> _loci_names;

	void _writeLociNames() {
		_lociNamesFile.numCols(_loci_names.size());
		for (auto &it : _loci_names) _lociNamesFile << it;
		_lociNamesFile.endln();
	};

	void _initOutputFiles() override {
		_lfmmFile.open(_outName + ".lfmm");
		_lociNamesFile.open(_outName + ".lfmm.kept_loci");
	};

	void _closeOutputFiles() {
		_lfmmFile.close();
		_lociNamesFile.close();
	};

	void _writeLFMM() {
		size_t numLoci = _genotypes.dimensions()[0];
		for (size_t i = 0; i < _samples.numSamples(); i++) {
			for (size_t l = 0; l < numLoci; l++) { _lfmmFile << (double)_genotypes(l, i); }
			_lfmmFile.endln();
		}
	}

	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override {
		// LFMM has individuals as rows and loci as columns (transposed)
		// -> we need to store these values first and then write
		_store(data);
		_loci_names.emplace_back(_reader.chrPos());
	}

	virtual void _store(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) = 0;

public:
	TVcfToLFMM() : TVcfConverter(){};
	~TVcfToLFMM() override = default;

	void run() {
		_genotypes.prepareFillData(10000, {_samples.numSamples()});
		_parseVCF();
		_genotypes.finalizeFillData();

		// write actual lfmm
		_lfmmFile.numCols(_genotypes.dimensions()[0]); // we only know now how many loci there are
		_writeLFMM();

		// write loci names
		_writeLociNames();

		_closeOutputFiles();
	}
};

//------------------------------------------
// TVcfToLFMMCalledGeno
//------------------------------------------
class TVcfToLFMMCalledGeno : public TVcfToLFMM<uint8_t> {
private:
	void _store(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;

public:
	TVcfToLFMMCalledGeno();
	~TVcfToLFMMCalledGeno() override = default;
};

//------------------------------------------
// TVcfToLFMMPostGeno
//------------------------------------------
class TVcfToLFMMPostGeno : public TVcfToLFMM<double> {
private:
	void _store(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;

public:
	TVcfToLFMMPostGeno();
	~TVcfToLFMMPostGeno() override = default;
};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
class TVcfToPosFile : public TVcfConverter {
	// posFile is needed as input for STITCH
	// format:
	//   - tab-separated
	//   - no header
	//   - 4 columns: col 1 = chromosome, col 2 = physical position (sorted from smallest to largest), col 3 =
	//   reference base, col 4 = alternate base
	//   - one row per SNP
	// Idea: one posfile per chromosome -> use option keepChromosomes

private:
	coretools::TOutputFile _posFile;

	void _initOutputFiles() override;
	void _writeHeader() override;
	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;

public:
	~TVcfToPosFile() override = default;
	void run();
};

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------
class TVcfToGenotypeTruthSetFile : public TVcfConverter {
	// produces genotype truth set (genfile) for STITCH and bed files for samples
	// idea: first locus -> find 0-n samples that have depth higher than minDepth
	//                   -> write position of this locus into bed-files for these individuals
	//                   -> write genotypes of these individuals to genfile; write genotypes of all other individuals
	//                   as NA to genfile
	//       next locus  -> is distance to previous locus more than x basepairs?
	//                   -> If yes: find 0-n samples that have depth higher than minDepth
	//                              -> write genotypes of these individuals to genfile; write genotypes of all other
	//                              individuals as NA to genfile
	//                   -> Else: write genotypes of all individuals as NA to genfile
	// format:
	//   - produces a BED file (3 columns, col 1 = chromosome, col 2 = start (1-based), col 3 = stop) for each
	//   individual
	//   - produces a genfile
	//      - tab-separated
	//      - header = sample names from vcf
	//      - one row per SNP
	//      - genotypes encoded as 0,1,2 or NA

private:
	std::vector<genometools::TBed> _bedFiles;
	coretools::TOutputFile _genFile;

	bool _first                        = true;
	size_t _minDistanceToPreviousLocus = 0;
	size_t _positionPreviousLocus      = 0;
	size_t _numSamplesPerLocus         = 0;
	std::string _curChr;

	void _writeHeader() override;
	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) override;
	void _filterIndividuals(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data);
	void _mapIndividualsToDepth(std::vector<size_t> &samplesToKeep);
	void _filterIndividualsWithHighestDepth(
	    std::vector<size_t> &samplesToKeep,
	    const std::map<double, std::vector<size_t>, std::greater<>> &depthVsSampleIndexMap) const;
	void _writeToGenFile(const std::vector<size_t> &samplesToKeep);
	void _storeInBedFile(const std::vector<size_t> &samplesToKeep);
	void _initOutputFiles() override;
	void _resetDistance();

public:
	TVcfToGenotypeTruthSetFile();
	~TVcfToGenotypeTruthSetFile() override = default;
	void run();
};

//--------------------------------------
// Tasks
//--------------------------------------
class TTask_VcfConverter : public coretools::TTask {
public:
	TTask_VcfConverter() { _explanation = "Converting a VCF file to other formats"; };

	void run() override {
		using namespace coretools::instances;
		std::string format = parameters().getParameter<std::string>("format");

		if (format == "beagle") {
			logfile().startIndent("Converting a VCF to Beagle format, old version (parameter 'format'):");
			TVcfToBeagle VcfToBeagle;
			VcfToBeagle.run();
		} else if (format == "beagleNew") {
			logfile().startIndent("Converting a VCF to Beagle format, new version (parameter 'format'):");
			TVcfBeagleNew VcfToBeagle;
			VcfToBeagle.run();
		} else if (format == "geno") {
			logfile().startIndent("Converting a VCF to geno format (parameter 'format'):");
			TVcfToGeno vcfToGeno;
			vcfToGeno.run();
		} else if (format == "LFMM") {
			logfile().startIndent("Converting a VCF to LFMM format (parameter 'format'):");

			// posterior or call?
			std::string genoType = parameters().getParameterWithDefault<std::string>("genotypes", "call");
			if (genoType == "posterior") {
				TVcfToLFMMPostGeno vcfToLFMMPostGeno;
				vcfToLFMMPostGeno.run();
			} else if (genoType == "call") {
				TVcfToLFMMCalledGeno vcfToLFMMCalledGeno;
				vcfToLFMMCalledGeno.run();
			} else {
				UERROR("Unknown genotype method '", genoType, "'! Use either 'call' or 'posterior'");
			}
		} else if (format == "posfile") {
			logfile().startIndent("Converting a VCF file to posfile format used by STITCH (parameter 'format'):");
			TVcfToPosFile VcfToPosFile;
			VcfToPosFile.run();
		} else if (format == "genfile") {
			logfile().startIndent(
			    "Converting a VCF file to genotype truth set format used by STITCH (parameter 'format'):");
			TVcfToGenotypeTruthSetFile VcfToGenotypeTruthSetFile;
			VcfToGenotypeTruthSetFile.run();
		} else {
			UERROR("Unknown format '", format,
				   "'! Use either 'beagle', 'geno', 'LFMM', 'posfile' or 'genfile'.");
		}
		logfile().endIndent();
	};
};

} // namespace VCF

#endif // ATLAS_TVCFCONVERTER_H
