#include "TGLFPrinter.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringManipulations.h"
#include "coretools/Strings/toString.h"
#include "genometools/GLF/TGLFReader.h"
#include "genometools/GLF/TSingleGLFTraverser.h"
#include "genometools/Genotypes/Genotype.h"

namespace PopulationTools {

using coretools::instances::logfile;

void TGLFPrinter::run() {
	using coretools::str::readBeforeLast;
	using coretools::str::toString;
	using genometools::Base;
	using genometools::Genotype;

	genometools::TBed bed;
	if (!_bedName.empty()) {
		bed.parse(_bedName);
		logfile().list("Limiting to sites marked in bed-file '", _bedName, "'. (parameter 'bed').");
	} else {
		logfile().list("Printing all sites. (Limit with 'bed').");
	}
	genometools::TSingleGLFTraverser traverser(_glfName, bed);

	const auto outName = toString(readBeforeLast(_glfName, ".glf"),  ".txt");
	logfile().list("Writing GLF-file to '", outName, "'.");

	coretools::TOutputFile out(outName);
	out.writeNoDelim("GLF version: ", traverser.glfVersion()).endln();

	for (; !traverser.endOfChrs(); traverser.nextChr()) {
		const auto& curChr = traverser.curChr();
		out.writeNoDelim("CHROMOSOME: '", curChr.name(), "', length = ", curChr.length(),
						 ", ploidy = ", curChr.ploidy()).endln();
		out.write("chr", "pos", "depth", "RMSMappingQual");
		if (curChr.isHaploid()) {
			for (auto b = Base::min; b < Base::max; ++b) out.writeNoDelim("LL_", b).writeDelim();
		} else {
			for (auto g = Genotype::min; g < Genotype::max; ++g) out.writeNoDelim("LL_", g).writeDelim();
		}
		out.endln();

		for (; !traverser.endOfCurChr(); traverser.nextSite()) {
			const auto &e = traverser.site();
			out.writeln(curChr.name(), traverser.position().position() + 1, e.depth, e.RMSMappingQual, e.likelihoods);
			
		}
		out.writeln("---");
	}
}
} // namespace PopulationTools
