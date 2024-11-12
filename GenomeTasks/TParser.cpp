#include "TParser.h"
#include "TGenome.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {
using coretools::instances::logfile;
using coretools::instances::parameters;

TParser::TParser() {
	if (parameters().exists("trim3") || parameters().exists("trim5")) {
		_trim3 = parameters().get<int>("trim3", 0);
		if (_trim3 < 0) UERROR("trimming distance trim3 must be >= 0!");
		_trim5 = parameters().get<int>("trim5", 0);
		if (_trim5 < 0) UERROR("trimming distance trim5 must be >= 0!");
		if (_trim3 > 0 || _trim5 > 0) {
			logfile().list("Will trim first ", _trim3, " bases from the 3' end and the first ", _trim5,
						   " bases from the 5' end. (parameters 'trim3', 'trim5', respectively)");
		}
		_trimReads = true;
	} else {
		_trim3     = 0;
		_trim5     = 0;
		_trimReads = false;
	}
}

void TParser::openReference(bool required) {
	if (!_reference.isOpen()) {
		if (parameters().exists("fasta")) {
			std::string fastaFile = parameters().get<std::string>("fasta");
			logfile().list("Reading reference sequence from '" + fastaFile + "'. (parameter fasta)");
			_reference.open(fastaFile);
		} else {
			if (required) { UERROR("No reference provided! (Use parameter fasta to provide a reference)"); }
		}
	}
}

void TParser::fill(const TGenome& genome, BAM::TAlignment& alignment) const {
	// parse
	genome.bamFile().fill(alignment);
	alignment.parse(genome.errorModels().sequencingErrorModels());

	if (_trimReads) { alignment.trimRead(_trim3, _trim5); }
	alignment.filter(_qualityFilter); // always on
	if (_contextFilter) alignment.filter(_contextFilter);
	if (_reference.isOpen()) { alignment.addReference(_reference); }
}

} // namespace GenomeTasks
