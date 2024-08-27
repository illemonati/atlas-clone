/*
 * TPileup.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TPileup.h"

#include <algorithm>
#include <vector>

#include "coretools/Containers/TStrongArray.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"

#include "genometools/Genotypes/Base.h"

#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {

void addNames(std::vector<std::string> &vec) {
	using GT = genometools::Genotype;
	for (auto g = GT::min; g < GT::max; ++g) vec.push_back("P(D|" + genometools::toString(g) + ")");
}

bool parseField(std::set<std::string> &fields, const std::string &tag, const std::string &explanation) {
	if (fields.find(tag) != fields.end()) {
		logfile().list(explanation, " (", tag, ")");
		fields.erase(fields.find(tag));
		return true;
	}
	return false;
}

} // namespace impl

//---------------------------------
// TPileup
//---------------------------------
TPileup::TPileup() {
	_onlySummary = parameters().exists("onlySummaries");
	if (!_onlySummary) {
		// open output file
		const std::string filename = _genome.front().outputName() + "_pileup.txt.gz";
		logfile().list("Writing pileup to file '", filename, "'.");
		_out.open(filename);

		// parse output fields
		logfile().startIndent("Will print the following pileup fields (parameter 'fields'):");
		const auto tmp = parameters().get<std::vector<std::string>>(
		    "fields", {"depth", "bases", "qualities", "alleles", "mates", "strands", "likelihoods", "hml"});
		std::set<std::string> fields(tmp.begin(), tmp.end());

		_printSettings.set<Print::Depth>(impl::parseField(fields, "depth", "Sequencing depth"));
		_printSettings.set<Print::Bases>(impl::parseField(fields, "bases", "Pileup bases"));
		_printSettings.set<Print::SampleBases>(impl::parseField(fields, "sampleBases", "Sample bases per bam-file"));
		_printSettings.set<Print::Qualities>(impl::parseField(fields, "qualities", "Pileup qualities"));
		_printSettings.set<Print::Alleles>(impl::parseField(fields, "alleles", "Allele counts"));
		_printSettings.set<Print::Mates>(impl::parseField(fields, "mates", "Mate information"));
		_printSettings.set<Print::Strand>(impl::parseField(fields, "strands", "Strand information"));
		_printSettings.set<Print::Likelihoods>(impl::parseField(fields, "likelihoods", "Genotype likelihoods"));
		_printSettings.set<Print::HML>(impl::parseField(fields, "hml", "Heterozygot is most likely"));
		logfile().endIndent();

		// check if unknown fields were given
		if (!fields.empty()) {
			if (fields.size() == 1) {
				UERROR("Unknown field '", *fields.begin(),
				       "'! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates' and 'strands'.");
			} else {
				std::string f;
				for (auto i : fields) { f += '\'' + i + "', "; }
				UERROR("Unknown fields: ", f.substr(0, f.size() - 2),
				       "! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates',  'strands' and "
				       "'likelihoods'.");
			}
		}

		// compile header
		std::vector<std::string> header = {"chromosome", "position"};
		if (_windows.parser().reference()) { header.push_back("reference"); }
		if (_printSettings.get<Print::Depth>()) {
			header.push_back("depth");
			if (_windows.parser().reference()) {
				header.push_back("numRef");
				header.push_back("numNonRef");
			}
		}
		if (_printSettings.get<Print::Bases>()) { header.push_back("bases"); }
		if (_printSettings.get<Print::SampleBases>()) {
			header.push_back("sampledBases");
			header.push_back("numSampledA");
			header.push_back("numSampledC");
			header.push_back("numSampledG");
			header.push_back("numSampledT");
			header.push_back("totSampled");
		}
		if (_printSettings.get<Print::Qualities>()) { header.push_back("qualities"); }
		if (_printSettings.get<Print::Alleles>()) {
			header.push_back("numA");
			header.push_back("numC");
			header.push_back("numG");
			header.push_back("numT");
			header.push_back("numAlleles");
		}
		if (_printSettings.get<Print::Mates>()) {
			header.push_back("numFirstMate");
			header.push_back("numSecondMate");
		}
		if (_printSettings.get<Print::Strand>()) {
			header.push_back("numForwardStrand");
			header.push_back("numReverseStrand");
		}
		if (_printSettings.get<Print::Likelihoods>()) {
			using GT = genometools::Genotype;
			for (auto g = GT::min; g < GT::max; ++g) header.push_back("P(D|" + genometools::toString(g) + ")");
			header.push_back("mostLikelyGenotype");
		}
		if (_printSettings.get<Print::HML>()) {
			header.push_back("hetMostLikely");
		}

		_out.writeHeader(header);

		// print all sites, also those without data?
		if (parameters().exists("printAll")) {
			_printSettings.set<Print::OnlySitesWithData>(false);
			logfile().list(
			    "Will print all sites that pass filters, including those without data. (parameter 'printAll')");
		} else {
			_printSettings.set<Print::OnlySitesWithData>(true);
			logfile().list("Will print only sites with data. (use 'printAll' to print all)");
		}
	}

	if (parameters().exists("histograms") || _onlySummary) {
		logfile().startIndent("Will print the following histograms (parameter 'histograms'):");
		const auto tmp = parameters().get<std::vector<std::string>>("histograms", {"depth", "qualities", "contexts"});
		std::set<std::string> histograms(tmp.begin(), tmp.end());

		if (histograms.empty()) {
			histograms.emplace("depth");
			histograms.emplace("qualities");
			histograms.emplace("contexts");
		}

		_histSettings.set<Hist::Depths>(impl::parseField(histograms, "depth", "Sequencing depth"));
		_histSettings.set<Hist::Quality>(impl::parseField(histograms, "qualities", "Base qualities"));
		_histSettings.set<Hist::Contexts>(impl::parseField(histograms, "contexts", "Base contexts"));
		_histSettings.set<Hist::AllelicDepth>(impl::parseField(histograms, "allelicDepth", "Allelic depth"));
		logfile().endIndent();

		// check if unknown fields were given
		if (!histograms.empty()) {
			if (histograms.size() == 1) {
				UERROR("Unknown histogram '", *histograms.begin(),
				       "'! Valid histograms are 'depth', 'qualities', 'allelicDepth' and 'contexts'.");
			} else {
				std::string f;
				for (auto i : histograms) { f += '\'' + i + "', "; }
				UERROR("Unknown histograms: ", f.substr(0, f.size() - 2),
				       "! Valid histograms are 'depth', 'qualities', 'allelicDepth' and 'contexts'.");
			}
		}
		if (_histSettings.get<Hist::Depths>()) {
			_outDepthHistogram.open(_genome.front().outputName() + "_depthPerWindow.txt.gz", {"window", "depth"});
			_outDepthPerChromosome.open(_genome.front().outputName() + "_depthPerChromosome.txt.gz", {"chromosome", "depth"});
		}

		if (_histSettings.get<Hist::AllelicDepth>()) {
			logfile().list("Will assemble allelic depth up to a max depth of ", _windows.uptoDepth(), ". (parameter 'readUpToDepth')");
			if (_windows.uptoDepth()> 100) {
				logfile().warning("Allocating count table for a max depth of ",_windows.uptoDepth(), " uses a lot of memory! Use argument readUpToDepth to limit.");
			}

			_counts.resize(_windows.uptoDepth());

			if (parameters().exists("includeZero")) {
				_writeEmpty = true;
				logfile().list("Will write full table, including cells with zero counts. (parameter 'includeZero')");
			} else {
				_writeEmpty = false;
				logfile().list("Will only print cells with non-zero counts. (use 'includeZero' to print all cells)");
			}
		}
	} else {
		logfile().list("Will not output histograms (use 'histograms' to do so).");
	}
}

void TPileup::_handleWindow(GenotypeLikelihoods::TWindow& window) {
	using genometools::Base;
	logfile().list("Writing pileup ...");

	if (_histSettings.get<Hist::AllelicDepth>()) { logfile().list("Adding sites to allelic depth table ..."); }
	for (size_t pos = 0; pos < window.size(); ++pos) {
		if (window.isMasked(pos)) continue;
		const auto &site = window[pos];

		// write histograms
		if (_histSettings.get<Hist::Depths>()) {
			_depthPerSite.add(site.depth());
			_depthPerSitePerChromosome.add(site.depth());
		}

		if (_histSettings.get<Hist::AllelicDepth>()) {
			const auto alleleCounts = site.countAlleles();
			_counts.addSite(alleleCounts);
		}

		if (_histSettings.get<Hist::Quality>()) {
			for (auto &b : site) {
				if (b.base != genometools::Base::N) { _qualDist.add(b.readGroupID, b.recalQuality.get()); }
			}
		}

		if (_histSettings.get<Hist::Contexts>()) {
			for (auto &b : site) { _contextDist.add(b.recalQuality.get(), coretools::index(b.context())); }
		}

		if ((_printSettings.get<Print::OnlySitesWithData>() && site.empty()) || _onlySummary) continue;
		_out.write(window.chrName(),
		           window.positionOnChr(pos) + 1); // positions are zero-based internally

		if (_windows.parser().reference()) { _out.write(site.refBase); }
		if (_printSettings.get<Print::Depth>()) {
			_out.write(site.depth());
			if (_windows.parser().reference()) { _out.write(site.refDepth(), site.depth() - site.refDepth()); }
		}
		if (_printSettings.get<Print::Bases>()) { _out.write(site.getBases()); }
		if (_printSettings.get<Print::SampleBases>()) {
			const auto sBases = site.sampleBases();
			if (sBases.empty()) {
				_out.write("-", 0, 0, 0, 0, 0);
			} else {
				coretools::TStrongArray<size_t, Base> counts{};
				size_t tot = 0;
				for (auto b : sBases) {
					_out.writeNoDelim(b);
					++counts[b];
					++tot;
				}
				_out.writeDelim();
				_out.write(counts, tot);
			}
		}
		if (_printSettings.get<Print::Qualities>()) { _out.write(site.getQualities()); }
		if (_printSettings.get<Print::Alleles>()) {
			const auto alleleCounts = site.countAlleles();
			_out.write(alleleCounts[Base::A], alleleCounts[Base::C], alleleCounts[Base::G], alleleCounts[Base::T]);
			_out.write((int)coretools::numNonZero(alleleCounts));
		}
		if (_printSettings.get<Print::Mates>()) {
			const auto mateCounts = site.countMates();
			_out.write(mateCounts);
		}
		if (_printSettings.get<Print::Strand>()) {
			const auto strandCounts = site.countFwdRev();
			_out.write(strandCounts);
		}
		if (_printSettings.get<Print::Likelihoods>() || _printSettings.get<Print::HML>()) {
			const auto genoLik = _genome.front().errorModels().calculateGenotypeLikelihoods(site);
			auto g = genometools::Genotype(std::max_element(genoLik.begin(), genoLik.end()) - genoLik.begin());
			if (_printSettings.get<Print::Likelihoods>()) { _out.write(genoLik, genometools::toString(g)); }
			if (_printSettings.get<Print::HML>()) { _out.write(genometools::isHeterozygous(g)); }
		}
		_out.endln();
	}

	// write depth per window
	// also write depth per chromosome if this window is the last window of a chromosome
	if (_histSettings.get<Hist::Depths>()) {
		logfile().list("Writing sequencing depth estimates to file ...");
		_outDepthHistogram
		    .writeNoDelim(window.chrName(), ':', window.from().position() + 1, '-', window.to().position())
		    .writeDelim();
		_outDepthHistogram.writeln(window.depth());
		logfile().done();
	}
}

void TPileup::_endChromosome(const genometools::TChromosome &Chr) {
	if (_histSettings.get<Hist::Depths>()) {
		_outDepthPerChromosome.writeln(Chr.name(), _depthPerSitePerChromosome.mean());
		_depthPerSitePerChromosome.clear();
	}
}

void TPileup::run() {
	_traverseBAMWindows();

	if (_histSettings.get<Hist::Depths>()) {
		// write distribution
		logfile().list("Writing depth per site distribution to file '", _genome.front().outputName(), "_depthPerSiteHistogram.txt.gz'");
		logfile().list("Writing average depth per window to file '", _genome.front().outputName(), "_depthPerWindow.txt.gz'");
		logfile().list("Writing average depth per chromosome to file '", _genome.front().outputName(), "_depthPerChromosome.txt.gz'");
		_depthPerSite.write(_genome.front().outputName() + "_depthPerSiteHistogram.txt.gz", "depth");
	}

	if (_histSettings.get<Hist::Quality>()) {
		// print distribution
		std::string outputFileName = _genome.front().outputName() + "_qualHistogram.txt.gz";
		logfile().list("Writing quality distribution to '", outputFileName, "'.");
		coretools::TOutputFile out(outputFileName, {"readGroup", "quality", "counts"});

		// get read group names
		std::vector<std::string> readGroupNames;
		_genome.front().bamFile().readGroups().fillVectorWithNames(readGroupNames);
		// write combined
		_qualDist.resize(readGroupNames.size()); // make sure it has the right size
		_qualDist.writeCombined(out, "allReadGroups");
		_qualDist.write(out, readGroupNames);
	}

	if (_histSettings.get<Hist::Contexts>()) {
		// write counts
		std::string outputFileName = _genome.front().outputName() + "_contextInformation.txt.gz";
		logfile().list("Writing context information to file '", outputFileName, "'.");

		std::vector<std::string> contextLabels;

		for (genometools::BaseContext c = genometools::BaseContext::min; c < genometools::BaseContext::max; ++c) {
			contextLabels.push_back(genometools::toString(c));
		}

		_contextDist.resize(contextLabels.size()); // make sure it has the right size
		_contextDist.writeAsMatrix(outputFileName, "quality", contextLabels);
	}

	if (_histSettings.get<Hist::AllelicDepth>()) {
		// write to file
		std::string outputFileName = _genome.front().outputName() + "_allelicDepth.txt.gz";
		logfile().list("Writing allelic depth table to '", outputFileName, "' ...");
		_counts.write(outputFileName, _writeEmpty);
		logfile().done();
	}
}

} // namespace GenomeTasks
