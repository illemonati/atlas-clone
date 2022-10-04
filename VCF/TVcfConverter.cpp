//
// Created by caduffm on 2/24/20.
//

#include "TVcfConverter.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <memory>
#include <ostream>
#include <string>
#include <string_view>

#include <fmt/os.h>

#include "genometools/GenotypeTypes.h"
#include "genometools/BED/TBed.h"
#include "genometools/VCF/TPopulationLikelihoodLocus.h"
#include "coretools/Files/gzstream.h"
#include "coretools/Types/probability.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Strings/fromString.h"

namespace VCF {

using namespace coretools::instances;

//------------------------------------------
// TVcfConverter
//------------------------------------------

TVcfConverter::TVcfConverter() {
	_openVCF();
	_readSamples();

	// read output name
	std::string tmp = coretools::str::readBeforeLast(_vcfName, ".vcf");
	_outName        = parameters().getParameterWithDefault<std::string>("out", tmp);
	logfile().list("Writing output files with prefix '" + _outName + "'. (specify with 'out')");
}

void TVcfConverter::_openVCF() {
	_vcfName = parameters().getParameterFilename("vcf");
	_reader.initialize(false);
	_reader.openVCF(_vcfName);
}

void TVcfConverter::_readSamples() {
	if (parameters().parameterExists("samples")) {
		_samples.readSamples(parameters().getParameter<std::string>("samples"));
	}

	// match samples
	if (_samples.hasSamples()) {
		_samples.fillVCFOrder(_reader.getSampleVCFNames());
	} else {
		_samples.readSamplesFromVCFNames(_reader.getSampleVCFNames());
	}
}

void TVcfConverter::_parseVCF() {
	// write header
	_initOutputFiles();
	_writeHeader();

	// initialize storage
	genometools::TPopulationLikehoodLocus<TSampleLikelihoods> data(_samples.numSamples());

	// run through VCF file
	logfile().startIndent("Parsing VCF file...");
	while (_reader.readDataFromVCF(data, _samples)) { _write(data); }
	logfile().endIndent();

	// end of vcf file reached
	_reader.concludeFilters();
}

void TVcfConverter::_writeHeader() {
	// can stay empty if there is no header
}

//------------------------------------------
// TVcfToBeagle
//------------------------------------------

TVcfToBeagle::TVcfToBeagle() : TVcfConverter() {}

void TVcfToBeagle::_initOutputFiles() { _beagleFile.open(_outName + ".beagle.gz"); }

void TVcfToBeagle::_writeHeader() {
	std::vector<std::string> header = {"marker", "allele1", "allele2"};
	header.reserve(3 + 3 * _samples.numSamples()); // 3 entries per sample (3 genotypes)
	for (size_t s = 0; s < _samples.numSamples(); s++) {
		std::string sampleName = _samples.sampleName(s);
		for (size_t i = 0; i < 3; i++) { header.push_back(sampleName); }
	}
	_beagleFile.writeHeader(header);
}

void TVcfToBeagle::_writeHaploid(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data, size_t s) {
	using coretools::Probability;
	using BG      = genometools::BiallelicGenotype;
	_foundHaploid = true;

	// FastNGSadmix: use 3 columns and set 3rd gtl for haploids to zero
	// (http://www.popgen.dk/software/index.php/FastNGSadmix)
	// Not sure how other ANGSD tools handle haploid - nothing mentioned in documentation!
	if (data[s].isMissing()) {
		_beagleFile << 0.5 << 0.5 << 0; // use 0.5 instead of 1 as ANGSD requires GTL to sum to one
	} else {
		std::array<double, 2> GTL = {(Probability)data[s][BG::haploidFirst], (Probability)data[s][BG::haploidSecond]};
		coretools::normalize(GTL);
		_beagleFile << GTL[0] << GTL[1] << 0.0;
	}
}

void TVcfToBeagle::_writeDiploid(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data, size_t s) {
	using coretools::Probability;
	using BG = genometools::BiallelicGenotype;

	// ANGSD: 3rd gtl for haploids must be zero
	if (data[s].isMissing()) {
		constexpr static double m = 1. / 3.;
		_beagleFile << m << m << m; // use 1/3 instead of 1 as ANGSD requires GTL to sum to one
	} else {
		std::array<double, 3> GTL = {(Probability)data[s][BG::homoFirst], (Probability)data[s][BG::het],
		                             (Probability)data[s][BG::homoSecond]};
		coretools::normalize(GTL);
		_beagleFile << GTL[0] << GTL[1] << GTL[2];
	}
}

void TVcfToBeagle::_write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) {
	_beagleFile << _reader.chr() + "_" + coretools::str::toString(_reader.position());
	_beagleFile << _reader.refAllele() << _reader.altAllele();

	// write line
	for (size_t s = 0; s < _samples.numSamples(); s++) {
		if (data[s].isHaploid()) {
			_writeHaploid(data, s);
		} else {
			_writeDiploid(data, s);
		}
	}

	_beagleFile.endln();
}

void TVcfToBeagle::run() {
	_parseVCF();

	if (_foundHaploid) {
		// Note: fastNGSAdmix allows for haploid sites (http://www.popgen.dk/software/index.php/FastNGSadmix),
		// but flag -haploid must be used.
		logfile().warning("Detected haploid sites in VCF: only filled two first GTL columns and set 3rd to zero. Use "
		                  "downstream tools at own responsibility: e.g. fastNGSAdmix requires flag -haploid.");
	}

	_beagleFile.close();
}

class TLineReader {
	std::istream *_stream; // non-owning
	std::string _line;     // current line
public:
	// see https://www.informit.com/articles/printerfriendly/1407357 for interface discussion

	// Not taking ownership of stream
	TLineReader(std::istream & Stream): _stream(&Stream) {
		// Read first line
		popFront();
	}

	bool empty() const noexcept {
		assert(_stream);
		return _stream->eof();
	}

	// IMPORTANT: as soon as popFront is called, this string_view will be invalidated
	std::string_view front() const noexcept {
		assert(!empty());
		return std::string_view{_line}; //does not make a copy
	};
	void popFront() {
		assert(_stream);
		assert(!empty());
		std::getline(*_stream, _line);
	};
};

template<typename Delim = char>
class TSplitter {
	static_assert(std::is_same_v<Delim, char> || std::is_same_v<Delim, std::string> || std::is_same_v<Delim, std::string_view>);
	std::string_view _sv;
	Delim _delim;
	size_t _count;

public:
	TSplitter(std::string_view Sv, Delim delim) : _sv(Sv), _delim(delim), _count(_sv.find(_delim)) {}

	bool empty() const noexcept { return _sv.empty(); }

	std::string_view front() const noexcept {
		assert(!empty());
		return _sv.substr(0, _count);
	}

	void popFront() noexcept {
		assert(!empty());
		if (_count == std::string_view::npos) {
			_sv.remove_prefix(_sv.size());
		} else {
			if constexpr (std::is_same_v<Delim, char>) _sv.remove_prefix(_count + 1);
			else _sv.remove_prefix(_count + _delim.size());

			_count = _sv.find(_delim); // will be npos for last element
		}
	}
};

template<typename Range>
void skip(Range& range, size_t nGaps = 1) {
	for (size_t _ = 0; _ < nGaps; ++_) range.popFront();
}

void TVcfBeagleNew::run() {
	const auto inName = parameters().getParameterFilename("vcf");
	const auto outName =
		parameters().getParameterWithDefault<std::string>("out", coretools::str::readBeforeLast(inName, ".vcf")) +
		".beagle.gz";

	std::unique_ptr<std::istream> istream;
	if (coretools::str::readAfterLast(inName, '.') == "gz")
		istream = std::make_unique<gz::igzstream>(inName.c_str());
	else
		istream = std::make_unique<std::ifstream>(inName);

	TLineReader lineReader(*istream);

	// Read info lines
	bool hasGL = false;
	bool hasPL = false;

	for (; lineReader.front().substr(0, 2) == "##"; lineReader.popFront()) {
		const auto line = lineReader.front();
		constexpr std::string_view glString = "##FORMAT=<ID=GL";
		constexpr std::string_view plString = "##FORMAT=<ID=PL";
		if (line.substr(0, glString.size()) == glString) hasGL = true;
		if (line.substr(0, plString.size()) == plString) hasPL = true;
	}

	if (!hasGL && !hasPL) UERROR("vcf file needs field GL or PL");

	// Print header
	if (lineReader.front().substr(0, 6) != "#CHROM") UERROR("vcf file needs header");

	// open gzFile
	auto ostream = gz::ogzstream(outName.c_str());

	fmt::memory_buffer buf;
	auto bIt = std::back_inserter(buf);

	fmt::format_to(bIt, "marker\tallele1\tallele2");

	TSplitter header{lineReader.front(), '\t'};
	skip(header, 9);

	for (; !header.empty(); header.popFront()) {
		const auto s = header.front();
		for (size_t _ = 0; _ < 3; ++_) {
			fmt::format_to(bIt, "\t{}", s);
		}
	}
	fmt::format_to(bIt, "\n");
	ostream.write(buf.data(), buf.size());

	// Lines
	lineReader.popFront();
	TSplitter line1{lineReader.front(), '\t'};
	skip(line1, 8);

	// find GL
	TSplitter format{line1.front(), ':'};
	int nGL = -1;
	for(; !format.empty(); format.popFront()) {
		++nGL;
		if (format.front() == "GL") break;
	}
	if (format.empty()) UERROR("FORMAT string neets GL");

	for(; !lineReader.empty(); lineReader.popFront()) {
		buf.clear();
		auto bIt = std::back_inserter(buf);
		TSplitter line{lineReader.front(), '\t'};
		fmt::format_to(bIt, line.front());

		line.popFront();
		fmt::format_to(bIt, "_{}", line.front()); // POS

		line.popFront(); // skip ID
		line.popFront();
		fmt::format_to(bIt, "\t{}", line.front()); // REF

		line.popFront();
		fmt::format_to(bIt, "\t{}", line.front()); // ALT

		line.popFront(); // skip QUAL
		line.popFront(); // skip FILTER
		line.popFront(); // skip INFO
		line.popFront(); // skip FORMAT

		line.popFront(); // First sample
		for(; !line.empty(); line.popFront()) {
			TSplitter sample{line.front(), ':'};
			skip(sample, nGL); // go to GL field
			if (sample.front()[0] == '.') {
				fmt::format_to(bIt, "\t0.333333\t0.333333\t0.333333"); 
				continue;
			}
			TSplitter gls_sv{sample.front(), ','};
			std::array<double, 3> gls;
			double tot = 1.;

			for (size_t i = 0; i < 3; ++i) { // haploid left as an exercice
				const auto sv = gls_sv.front();
				if (sv[0] == '0') {
					// sample.front(): 0,-3.902,-17.902
					gls[i] = 1.;
				} else {
					coretools::str::fromString(sv, gls[i]);
					gls[i] = exp10(gls[i]);
					tot += gls[i];
				}
				gls_sv.popFront();
			}

			for (auto &gl : gls) { gl /= tot; }

			fmt::format_to(bIt, "\t{:.4}\t{:.4}\t{:.4}", gls[0], gls[1], gls[2]);
		}
		fmt::format_to(bIt, "\n");
		ostream.write(buf.data(), buf.size());
	}
}

//------------------------------------------
// TVcfToGeno
//------------------------------------------

TVcfToGeno::TVcfToGeno() : TVcfConverter() {}

void TVcfToGeno::_initOutputFiles() {
	_genoFile.open(_outName + ".geno");
	_lociNamesFile.open(_outName + ".geno.kept_loci");
}

void TVcfToGeno::_closeOutputFiles() {
	_genoFile.close();
	_lociNamesFile.endln();
	_lociNamesFile.close();
}

void TVcfToGeno::_write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) {
	_lociNamesFile << _reader.chrPos();

	// write line
	std::string line;
	for (size_t s = 0; s < _samples.numSamples(); s++) {
		if (data[s].isMissing()) {
			line += "9";
		} else if (data[s].isHaploid()) {
			line += coretools::str::toString(1 - altAlleleCounts(_reader.biallelicGenotype(_samples, s)));
		} else {
			line += coretools::str::toString(2 - altAlleleCounts(_reader.biallelicGenotype(_samples, s)));
		}
	}
	_genoFile.writeln(line);
}

void TVcfToGeno::run() {
	_parseVCF();
	_closeOutputFiles();
}

//------------------------------------------
// TVcfToLFMMCalledGeno
//------------------------------------------

TVcfToLFMMCalledGeno::TVcfToLFMMCalledGeno() : TVcfToLFMM<uint8_t>() {
	logfile().list("Will store the called genotype for each locus.");
}

void TVcfToLFMMCalledGeno::_store(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &) {
	// store allele counts
	std::vector<genometools::BiallelicGenotype> geno = _reader.biallelicGenotypes(_samples);
	for (size_t i = 0; i < _samples.numSamples(); i++) {
		if (isMissing(geno[i])) {
			_genotypes.emplace_back(9); // re-code missing genotypes to LFMM format
		} else {
			_genotypes.emplace_back(altAlleleCounts(geno[i]));
		}
	}
}

//------------------------------------------
// TVcfToLFMMPostGeno
//------------------------------------------

TVcfToLFMMPostGeno::TVcfToLFMMPostGeno() : TVcfToLFMM<double>() {
	logfile().list("Will store the mean posterior genotype for each locus.");
}

void TVcfToLFMMPostGeno::_store(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) {
	for (size_t i = 0; i < _samples.numSamples(); i++) {
		if (data[i].isMissing()) {
			UERROR("Missing data at sample ", _samples.sampleName(i), " and locus ", _reader.chr(), ":",
			       _reader.position(), "! LFMM2 does not accept missing genotypes, impute your VCF file first.");
		}
		_genotypes.emplace_back(data[i].meanPosteriorGenotype());
	}
}

//------------------------------------------
// TVcfToPosFile
//------------------------------------------

void TVcfToPosFile::_writeHeader() {
	// no header
	_posFile.numCols(4);
}

void TVcfToPosFile::_initOutputFiles() { _posFile.open(_outName + ".pos"); }

void TVcfToPosFile::_write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &) {
	_posFile.writeln(_reader.chr(), _reader.position(), _reader.refAllele(), _reader.altAllele());
}

void TVcfToPosFile::run() {
	// read Vcf and write output
	_parseVCF();

	// clean up
	_posFile.close();
}

//------------------------------------------
// TVcfToGenotypeTruthSetFile
//------------------------------------------

TVcfToGenotypeTruthSetFile::TVcfToGenotypeTruthSetFile() {
	// read params
	_minDistanceToPreviousLocus = parameters().getParameterWithDefault<size_t>("minDistance", 100);
	logfile().list("Will keep loci that have a minimal distance of ", _minDistanceToPreviousLocus,
	               " to previous locus (parameter 'minDistance').");
	_resetDistance();
	_numSamplesPerLocus = parameters().getParameterWithDefault<size_t>("numSamples", 5);
	logfile().list("Will keep up to ", _numSamplesPerLocus, " individuals per locus (parameter 'numSamples').");
}

void TVcfToGenotypeTruthSetFile::_writeHeader() { _genFile.writeHeader(_samples.sampleNames()); }

void TVcfToGenotypeTruthSetFile::_initOutputFiles() {
	_genFile.open(_outName + ".gen");

	// initialize bed files (we know how many samples there are)
	_bedFiles.resize(_samples.numSamples());
}

void TVcfToGenotypeTruthSetFile::_mapIndividualsToDepth(std::vector<size_t> &samplesToKeep) {
	std::map<double, std::vector<size_t>, std::greater<>> depthVsSampleIndexMap;
	// use depth as key in map -> automatically sorted in descending order
	for (auto &s : samplesToKeep) {
		// does depth already exist in map?
		double depth = _reader.depth(_samples, s);
		auto it      = depthVsSampleIndexMap.find(depth);
		if (it == depthVsSampleIndexMap.end()) { // new depth
			depthVsSampleIndexMap[depth] = {s};
		} else { // depth already exists
			it->second.push_back(s);
		}
	}
	_filterIndividualsWithHighestDepth(samplesToKeep, depthVsSampleIndexMap);
}

void TVcfToGenotypeTruthSetFile::_filterIndividualsWithHighestDepth(
    std::vector<size_t> &samplesToKeep,
    const std::map<double, std::vector<size_t>, std::greater<>> &depthVsSampleIndexMap) const {
	// only keep top _numSamplesPerLocus
	samplesToKeep.clear();
	size_t c = 0;
	for (auto &it : depthVsSampleIndexMap) { // one depth
		for (auto &it2 : it.second) {        // one sample at this depth
			if (c >= _numSamplesPerLocus) return;
			samplesToKeep.push_back(it2);
			c++;
		}
	}
}

void TVcfToGenotypeTruthSetFile::_filterIndividuals(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) {
	std::vector<size_t> samplesToKeep;

	size_t distanceToPreviousLocus = _reader.position() - _positionPreviousLocus;
	if (distanceToPreviousLocus >= _minDistanceToPreviousLocus || _first) { // check if distance is big enough
		// idea: TPopulationLikelihoods will filter on minDepth and set all samples with < minDepth as missing
		// here, we check how many individuals have > minDepth; we rank them and only keep _numSamplesPerLocus of them
		for (size_t s = 0; s < _samples.numSamples(); ++s) {
			if (!data[s].isMissing()) { samplesToKeep.push_back(s); }
		}
		if (samplesToKeep.empty()) { // no individuals have > minDepth
			_writeToGenFile(samplesToKeep);
			// no need to write to bed file
		} else if (samplesToKeep.size() <= _numSamplesPerLocus) { // keep all of them
			_writeToGenFile(samplesToKeep);
			_storeInBedFile(samplesToKeep);
		} else { // rank according to depth
			_mapIndividualsToDepth(samplesToKeep);
			_writeToGenFile(samplesToKeep);
			_storeInBedFile(samplesToKeep);
		}
	} else {
		_writeToGenFile(samplesToKeep);
		// no need to write to bed file
	}
}

void TVcfToGenotypeTruthSetFile::_writeToGenFile(const std::vector<size_t> &samplesToKeep) {
	for (size_t s = 0; s < _samples.numSamples(); s++) {
		// should we write true genotype of sample?
		auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
		if (it != samplesToKeep.end()) {
			// sample found
			_genFile << toString(_reader.genotype(_samples, s));
		} else {
			_genFile << "NA";
		}
	}
	_genFile.endln();
}

void TVcfToGenotypeTruthSetFile::_storeInBedFile(const std::vector<size_t> &samplesToKeep) {
	for (size_t s = 0; s < _samples.numSamples(); s++) {
		// should we write to bed of sample?
		auto it = std::find(samplesToKeep.begin(), samplesToKeep.end(), s);
		if (it != samplesToKeep.end()) // sample found
			_bedFiles[s].add(_reader.chr(), _reader.position());
	}
	_positionPreviousLocus = _reader.position();
}

void TVcfToGenotypeTruthSetFile::_write(genometools::TPopulationLikehoodLocus<TSampleLikelihoods> &data) {
	if (_curChr != _reader.chr()) { // new chromosome
		_curChr = _reader.chr();
		_resetDistance();
	} else {
		_first = false;
	}
	_filterIndividuals(data);
}

void TVcfToGenotypeTruthSetFile::_resetDistance() {
	_positionPreviousLocus = 0; // for first locus on chromosome
	_first                 = true;
}

void TVcfToGenotypeTruthSetFile::run() {
	// read Vcf and write output
	_parseVCF();

	// write bed files (one per sample)
	for (uint32_t s = 0; s < _samples.numSamples(); s++) {
		// check if sample name contains / (would be interpreted as path)
		std::string sample_name = _samples.sampleName(s);
		if (coretools::str::stringContains(sample_name, '/'))
			sample_name = coretools::str::stringReplace("/", "_", sample_name);
		_bedFiles[s].write(_outName + "_" + sample_name + ".bed");
	}

	// clean up
	_genFile.close();
}

} // namespace VCF
