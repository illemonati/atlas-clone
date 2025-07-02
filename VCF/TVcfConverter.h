//
// Created by caduffm on 2/24/20.
//

#ifndef ATLAS_TVCFCONVERTER_H
#define ATLAS_TVCFCONVERTER_H

#include <map>
#include <string>
#include <vector>

#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TError.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/TSampleLikelihoods.h"
#include "genometools/VCF/TPopulationLikelihoods.h"

namespace genometools {
template<typename Type> class TPopulationLikehoodLocus;
}

namespace VCF {

using TSampleLikelihoods = genometools::TSampleLikelihoods<coretools::HPPhredInt>;

//------------------------------------------
// TVcfConverter
//------------------------------------------
class TVcfFileConverter {
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
	TVcfFileConverter();
	virtual ~TVcfFileConverter() = default;
};

//------------------------------------------
// TVcfToBeagle
//------------------------------------------
class TVcfToBeagle : protected TVcfFileConverter {
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
class TVcfToGeno : protected TVcfFileConverter {
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
// TVcfTranspose
//------------------------------------------

template<bool StoreCalledGenotypes> class TVcfTranspose : protected TVcfFileConverter {
	/*
	 * Base class for all file formats that require transposing
	 * i.e. rows are individuals and columns are loci (-> transposed compared to VCF format)
	 * this requires all genetic data to be stored in memory
	 * template parameter StoreCalledGenotypes:
	 * whether to store called genotypes or store mean posterior genotype
	 */
	using Type = std::conditional_t<StoreCalledGenotypes, uint8_t, double>;

protected:
	// store genetic data (called genotypes or mean posterior genotype)
	coretools::TMultiDimensionalStorage<Type, 2> _genotypes;

	// store chr:pos names
	std::vector<std::string> _loci_names;

	void _storeAlleleCounts() {
		std::vector<genometools::BiallelicGenotype> geno = _reader.biallelicGenotypes(_samples);
		for (size_t i = 0; i < _samples.numSamples(); i++) {
			isMissing(geno[i]) ? _genotypes.emplace_back(9) : _genotypes.emplace_back(altAlleleCounts(geno[i]));
		}
	}

	void _storeMeanPosteriorGenotype(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &Data) {
		// store mean posterior genotype
		for (size_t i = 0; i < _samples.numSamples(); i++) {
			coretools::user_assert(!Data[i].isMissing(), "Missing data at sample ", _samples.sampleName(i),
								   " and locus ", _reader.chr(), ":", _reader.position(), "!");
			_genotypes.emplace_back(Data[i].meanPosteriorGenotype());
		}
	}

	void _write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &Data) override {
		// store locus name
		_loci_names.emplace_back(_reader.chrPos());

		// store genetic data
		if constexpr (StoreCalledGenotypes) {
			_storeAlleleCounts();
		} else {
			_storeMeanPosteriorGenotype(Data);
		}
	}

	virtual void _writeTransposedFiles() = 0;

public:
	TVcfTranspose() : TVcfFileConverter(){};

	~TVcfTranspose() override = default;

	void run() {
		_genotypes.prepareFillData(10000, {_samples.numSamples()});
		_parseVCF();
		_genotypes.finalizeFillData();

		// write actual file
		_writeTransposedFiles();
	}
};

//------------------------------------------
// TVcfToLFMM
//------------------------------------------

template<bool StoreCalledGenotypes> class TVcfToLFMM : public TVcfTranspose<StoreCalledGenotypes> {
protected:
	using TVcfTranspose<StoreCalledGenotypes>::_loci_names;
	using TVcfTranspose<StoreCalledGenotypes>::_genotypes;
	using TVcfTranspose<StoreCalledGenotypes>::_outName;
	using TVcfTranspose<StoreCalledGenotypes>::_samples;
	using TVcfTranspose<StoreCalledGenotypes>::_reader;

	// files
	coretools::TOutputFile _lfmmFile;
	coretools::TOutputFile _lociNamesFile;

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

	void _writeTransposedFiles() override {
		// write actual lfmm
		_lfmmFile.numCols(_genotypes.dimensions()[0]); // we only know now how many loci there are
		_writeLFMM();

		// write loci names
		_writeLociNames();

		_closeOutputFiles();
	}

public:
	TVcfToLFMM() : TVcfTranspose<StoreCalledGenotypes>() {
		if constexpr (StoreCalledGenotypes) {
			coretools::instances::logfile().list("Will store the called genotype for each locus.");
		} else {
			coretools::instances::logfile().list("Will store the mean posterior genotype for each locus.");
		}
	};
	~TVcfToLFMM() override = default;
};

//------------------------------------------
// TVcfToSambada
//------------------------------------------
class TVcfToSambada : public TVcfTranspose<true> {
private:
	coretools::TOutputFile _sambadaFile;

	void _writeSambadaHeader();
	void _initOutputFiles() override;
	void _writeTransposedFiles() override;

public:
	TVcfToSambada();
	~TVcfToSambada() override = default;
};

//------------------------------------------
// TVcfToPosFile
//------------------------------------------
class TVcfToPosFile : public TVcfFileConverter {
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
class TVcfToGenotypeTruthSetFile : public TVcfFileConverter {
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

struct TVCFConverter {
	void run() {
		using coretools::instances::parameters;
		using coretools::instances::logfile;
		std::string format = parameters().get<std::string>("format");

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
			std::string genoType = parameters().get<std::string>("genotypes", "call");
			if (genoType == "posterior") {
				TVcfToLFMM<false> vcfToLFMMPostGeno;
				vcfToLFMMPostGeno.run();
			} else if (genoType == "call") {
				TVcfToLFMM<true> vcfToLFMMCalledGeno;
				vcfToLFMMCalledGeno.run();
			} else {
				throw coretools::TUserError("Unknown genotype method '", genoType, "'! Use either 'call' or 'posterior'");
			}
		} else if (format == "Sambada") {
			logfile().startIndent("Converting a VCF to Sambada format (parameter 'format'):");
			TVcfToSambada vcfToSambada;
			vcfToSambada.run();
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
			throw coretools::TUserError("Unknown format '", format, "'! Use either 'beagle', 'geno', 'LFMM', 'posfile' or 'genfile'.");
		}
		logfile().endIndent();
	};
};

} // namespace VCF

#endif // ATLAS_TVCFCONVERTER_H
