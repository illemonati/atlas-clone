#include "TOutputBamFile.h"

#include "TBamFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Math/TNumericRange.h"
#include "coretools/Types/probability.h"

namespace BAM {
using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::TNumericRange;
using coretools::PhredInt;

TQualityAdjusterForWriting::TQualityAdjusterForWriting() {
	if (parameters().exists("outQual")) {
		TNumericRange<char> qualRange;
		parameters().fill("outQual", qualRange);
		limitRange(qualRange);

		logfile().list("Will print qualities truncated to ", rangeString(), " (parameter 'outQual')");

		if (qualRange.max() > _maxBaseQuality) { logfile().warning("Truncated quality range to BAM limits!"); }

	} else {
		logfile().list(
			"Will use the full range of quality scores when writing alignments. (use 'outQual' to constrain).");
	}

	// quality binning
	if (parameters().exists("writeBinnedQualities")) {
		logfile().list("Will write Illumina-binned quality scores. (parameter 'writeBinnedQualities')");
		binQualitiesIllumina();
	} else {
		logfile().list("Will write raw quality scores. (use 'writeBinnedQualities' to bin)");
	}
};

void TQualityAdjusterForWriting::binQualitiesIllumina(){
	_binIllumina = true;
	_adjust      = true;
};

void TQualityAdjusterForWriting::limitRange(const TNumericRange<char> & Range){
	if(Range.minIncluded()){
		_minQual = Range.min();
	} else {
		_minQual = Range.min() + 1;
	}
	if(Range.maxIncluded()){
		_maxQual = Range.min();
	} else {
		_maxQual = Range.min() - 1;
	}
	_adjust = true;
};

std::string TQualityAdjusterForWriting::rangeString(){
	return "[" + coretools::str::toString(coretools::fromChar(_minQual)) + "," + coretools::str::toString(coretools::fromChar(_maxQual)) + "]";
};

char TQualityAdjusterForWriting::_adjustOneQuality(char Qual) const {
	if(Qual < _minQual){
		Qual = _minQual;
	} else if (Qual > _maxQual){
		Qual = _maxQual;
	}
	if(_binIllumina){
		Qual = makeIllumina(Qual);
	}

	return Qual;
};

void TQualityAdjusterForWriting::adjustQualities(std::string & qualities) const {
	if(_adjust){
		for(auto& q : qualities){
			q = _adjustOneQuality(q);
		}
	}
};

TOutputBamFile::TOutputBamFile(std::string Filename, const TBamFile &Original)
	: TOutputBamFile(std::move(Filename), Original.samHeader(), Original.chromosomes(), Original.readGroups()){};


TOutputBamFile::TOutputBamFile(std::string Filename, const TSamHeader &Header,
							   const genometools::TChromosomes &Chromosomes, const TReadGroups &ReadGroups,
							   const TQualityAdjusterForWriting &QualityAdjuster)
	: TOutputBamFile(std::move(Filename), Header, Chromosomes, ReadGroups) {
	_qualityAdjuster = QualityAdjuster;
}

TOutputBamFile::~TOutputBamFile(){
	// write all future alignments
	std::sort(_futureAlignments.begin(), _futureAlignments.end(), std::less<>());
	for (auto &a : _futureAlignments) { _writeAlignment(a); }
	_futureAlignments.clear();

	// close
	_bamWriter.Close();
	logfile().listFlush("Creating index of BAM file '" + _outputFilename + "' ...");

	// create index of BAM file
	BamTools::BamReader reader;
	if (!reader.Open(_outputFilename)) logfile().error("Failed to open BAM file '", _outputFilename, "' for indexing!");
	reader.CreateIndex(BamTools::BamIndex::STANDARD);

	// close BAM file
	reader.Close();
	logfile().done();
}

	TOutputBamFile::TOutputBamFile(std::string Filename, const TSamHeader & Header, const genometools::TChromosomes & Chromosomes, const TReadGroups & ReadGroups) :
		_outputFilename(std::move(Filename)), _readGroups(&ReadGroups){

	//construct new header /without chromosomes
	const auto header = Header.compileSamHeader(ReadGroups, Chromosomes);

	//fill bamtools chromosomes
	BamTools::RefVector ref;
	for(auto it = Chromosomes.cbegin(); it != Chromosomes.cend(); ++it){
		ref.emplace_back(it->name(), it->length());
	}

	//open file for writing
	if(!_bamWriter.Open(_outputFilename, header, ref)){
		UERROR("Failed to open BAM file '", _outputFilename, "'!");
	}
};

void TOutputBamFile::writeAlignment(BamTools::BamAlignment & alignment){
	//adjust qualities for printing
	_qualityAdjuster.adjustQualities(alignment.Qualities);

	// write alignment
	if(!_bamWriter.SaveAlignment(alignment))
		UERROR("Read '", alignment.Name, "' could not be written!");
};

void TOutputBamFile::_writeAlignment(const TAlignment & alignment){
	//create bamAlignment and then write
	BamTools::BamAlignment _tmpBamAlignment;

	_tmpBamAlignment.Name = alignment.name();
	_tmpBamAlignment.AlignmentFlag = alignment.flags();
	_tmpBamAlignment.RefID = alignment.refID();
	_tmpBamAlignment.Position = alignment.position();
	_tmpBamAlignment.InsertSize = alignment.insertSize();
	_tmpBamAlignment.MapQuality = alignment.mappingQuality().get();

	if(alignment.isPaired()){
		_tmpBamAlignment.MateRefID = alignment.mateRefID();
		_tmpBamAlignment.MatePosition = alignment.matePosition();
	}

	//CIGAR
	for(const auto& it : alignment.cigar()){
		_tmpBamAlignment.CigarData.emplace_back(it.type, it.length);
	}

	//add sequences and qualities
	_tmpBamAlignment.QueryBases = alignment.sequence();
	_tmpBamAlignment.Qualities = alignment.qualities();

	//add read group information
	_tmpBamAlignment.AddTag("RG", "Z", _readGroups->getName(alignment.readGroupId()));

	//and now write
	if(!_bamWriter.SaveAlignment(_tmpBamAlignment)){
		UERROR("Read '", _tmpBamAlignment.Name, "' could not be written!");
	}
};

void TOutputBamFile::writeAlignment(const TAlignment & alignment){
	//write alignments BEFORE alignment to write next
	std::sort(_futureAlignments.begin(), _futureAlignments.end(), std::less<>());
	std::reverse(_futureAlignments.begin(), _futureAlignments.end());
	while (!_futureAlignments.empty() && _futureAlignments.back().from() <= alignment.from()) {
		_writeAlignment(_futureAlignments.back());
		_futureAlignments.pop_back();
		}

	//write next alignment
	_writeAlignment(alignment);
};

void TOutputBamFile::writeAlignmentLater(const TAlignment & alignment) {
	_futureAlignments.push_back(alignment);
}
}
