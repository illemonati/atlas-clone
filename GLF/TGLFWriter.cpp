#include "TGLFWriter.h"
#include "TGLF.h"
#include "coretools/Files/TWriter.h"
#include "genometools/VCF/TVcfWriter.h"
#include <memory>


namespace GLF {

void TGLFWriter::_writeHeader(std::string_view Header) {
	_writer->write(version());

	_header = Header;
	if (_header.length() > 0) {
		const uint32_t labelLength = _header.size();
		_writer->write(labelLength);
		_writer->write(_header);
	} else {
		uint32_t zero32 = 0;
		_writer->write(zero32);
	}
};

void TGLFWriter::open(std::string_view Filename, const genometools::TChromosomes &, std::string_view Header){
	_writer.reset(coretools::makeWriter(Filename, "w"));

	if (!_writer->isOpen()) UERROR("Failed to open GLF file '", Filename, "' for writing!");
	_index.clear();

	// write header
	_writeHeader(Header);
};

void TGLFWriter::newChromosome(const genometools::TChromosome &chromosome) {
	// add chromosome info to index file
	_index.addChromosme(chromosome, _writer->tell());

	// write record type
	const uint8_t zero8 = 0;
	_writer->write(zero8);

	// write new chromosome: length of label, label, refId, length of ref sequence, ploidy
	const uint32_t labelLength = chromosome.name().size();
	_writer->write(labelLength);
	_writer->write(chromosome.name());
	_writer->write(chromosome.refID());
	_writer->write(chromosome.length());
	_writer->write(chromosome.ploidy()); // TODO: I get an "uninitialized variable" error with valgrind. Why?

	// set oldPos
	_oldPos = 0;
};

void TGLFWriter::writeSite(long pos, uint32_t depth, uint8_t RMS_mappingQual,
			   const GenotypeLikelihoods::TGenotypeLikelihoods &genotypeLikelihoods) {
	using genometools::Genotype;
	using coretools::Probability;
	using genometools::Ploidy;
	const uint8_t _recordType1 = 1 << 4;
	// record type
	// TODO: add reference?
	_writer->write(_recordType1);

	// offset
	const uint32_t offset = pos - _oldPos;
	_oldPos = pos;
	_writer->write(offset);

	TGLFLikelihoods glfValues; // tmp used for writing

	// calculate likelihoods in GLF format
	// Note: genotype likelihoods are given for the 10 diploid genotypes!!
	if (_index.chromosomes().back().isHaploid()) {
		using genometools::Base;
		const double maxLik = std::max({genotypeLikelihoods[Genotype::AA], genotypeLikelihoods[Genotype::CC],
										genotypeLikelihoods[Genotype::GG], genotypeLikelihoods[Genotype::TT]});

		// normalize and scale to uint16
		glfValues.type = Ploidy::haploid;
		glfValues[Base::A] = Probability(genotypeLikelihoods[Genotype::AA] / maxLik);
		glfValues[Base::C] = Probability(genotypeLikelihoods[Genotype::CC] / maxLik);
		glfValues[Base::G] = Probability(genotypeLikelihoods[Genotype::GG] / maxLik);
		glfValues[Base::T] = Probability(genotypeLikelihoods[Genotype::TT] / maxLik);
	} else {
		// ploidy is 2
		glfValues.type = Ploidy::diploid;
		const double maxLik = *std::max_element(genotypeLikelihoods.begin(), genotypeLikelihoods.end());

		// normalize and scale to genometools::HighPrecisionPhredIntProbability

		for (auto g = Genotype::min; g < Genotype::max; ++g) {
			glfValues[g] = Probability(genotypeLikelihoods[g]/maxLik);
		}
	}

	// write maxLL as uint16_t
	// uint16_t maxLL_int = converter.toGlfFormat(maxLL);
	// write(&maxLL_int, sizeof(uint16_t));

	// write depth as uint16_t
	depth = std::min(depth, uint32_t{65535});
	const uint16_t tmp = depth;
	_writer->write(tmp);

	// root mean square of mapping qualities
	_writer->write(RMS_mappingQual);

	// genotype likelihoods
	_writer->write(glfValues);
}
void TGLFWriter::close() {
	if(_writer->isOpen()){
		_index.writeChromosmes(_writer->name());
	}
	_writer = std::make_unique<coretools::TNoWriter>();
};
}
