/*
 * TPileup.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TPileup.h"

#include <algorithm>
#include <cmath>
#include <vector>


#include "coretools/Containers/TStrongArray.h"
#include "coretools/Containers/TView.h"
#include "coretools/Files/TOutputFile.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/algorithms.h"

#include "genometools/Genotypes/Base.h"

#include "TSequencedData.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks {

using coretools::instances::logfile;
using coretools::instances::parameters;
using coretools::TOutputFile;

namespace impl {

bool parseField(std::set<std::string> &fields, const std::string &tag, const std::string &explanation) {
	if (fields.find(tag) != fields.end()) {
		logfile().list(explanation, " (", tag, ")");
		fields.erase(fields.find(tag));
		return true;
	}
	return false;
}

template<typename Transitions>
void writeTransitions(const Transitions &transitions, std::string_view Chr, TOutputFile& _outTransitions,
					   TOutputFile& _outTransitionsRel, TOutputFile& _outTransitionsPsi, TOutputFile& _outTransitionsRho) {
	using BAM::End;
	using BAM::Mate;
	using BAM::Strand;
	using genometools::Base;
	using genometools::base2char;
	for (auto mate = Mate::min; mate < Mate::max; ++mate) {
		for (auto strand = Strand::min; strand < Strand::max; ++strand) {
			for (auto end = End::min; end < End::max; ++end) {
				coretools::TStrongArray<coretools::TStrongArray<size_t, genometools::Base>, genometools::Base> rho{};
				auto &tr = transitions[mate][strand][end];
				for (size_t i = 0; i < tr.size(); ++i) {
					_outTransitions.write(Chr, mate, strand, end, i);
					_outTransitionsRel.write(Chr, mate, strand, end, i);
					_outTransitionsPsi.write(Chr, mate, strand, end, i);
					_outTransitionsRho.write(Chr, mate, strand, end, i);
					coretools::TStrongArray<size_t, genometools::Base> tot{};
					for (auto ref = Base::min; ref < Base::max; ++ref) {
						// counts
						for (auto b = Base::min; b < Base::max; ++b) {
							_outTransitions.write(tr[i][ref][b]);
							tot[ref]       += tr[i][ref][b];
							rho[ref][b]    += tr[i][ref][b];
						}
						const auto totRho = tot[ref] - tr[i][ref][ref];

						for (auto b = Base::min; b < Base::max; ++b) {
							//Rel
							const auto rel = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
							_outTransitionsRel.write(fmt::format("{:.4f}", rel));

							//Rho
							if (ref == b) {
								_outTransitionsRho.write("  -  ");
							}
							else {
								const auto rho = totRho ? double(tr[i][ref][b])/totRho : 0;
								_outTransitionsRho.write(fmt::format("{:.4f}", rho));
							}
						}
					}
					// Psi
					for (auto ref = Base::min; ref < Base::max; ++ref) {
						for (auto b = Base::min; b < Base::max; ++b) {
							if (ref == b) {
								_outTransitionsPsi.write("0.000");
							} else {
								const auto fromTo = tot[ref] ? double(tr[i][ref][b]) / tot[ref] : 0;
								const auto toFrom = tot[b] ? double(tr[i][b][ref]) / tot[b] : 0;
								const auto Psi    = (fromTo - toFrom) / (1.0 - toFrom);
								_outTransitionsPsi.write(fmt::format("{:.4f}", Psi));
							}
						}
					}
					_outTransitions.endln();
					_outTransitionsRel.endln();
					_outTransitionsRho.endln();
					_outTransitionsPsi.endln();
				}
			}
		}
	}

	// Flush once per chromosome
	_outTransitions.flush();
	_outTransitionsRel.flush();
	_outTransitionsRho.flush();
	_outTransitionsPsi.flush();
}

} // namespace impl

//---------------------------------
// TPileup
//---------------------------------
TPileup::TPileup() {
	_onlyHistograms = parameters().exists("onlyHistograms");

	if (!_onlyHistograms) {
		// open output file
		const std::string filename = _siteTraverser.outputName() + "_pileup.txt.gz";
		logfile().list("Writing pileup to file '", filename, "'.");
		_out.open(filename);

		// parse output fields

		constexpr coretools::TStrongArray<std::string_view, Print> printNames{
			{"depth", "bases", "sampleBases", "qualities", "alleles", "mates", "strands", "likelihoods", "hml"}};
		constexpr coretools::TStrongArray<std::string_view, Print> printExpl{
			{"Sequencing depth", "Pileup bases", "Sample bases per bam-file", "Pileup qualities", "Allele counts",
			 "Mate information", "Strand information", "Genotype likelihoods", "Heterozygot is most likely"}};

		const auto fields = parameters().get<std::vector<std::string_view>>(
		    "fields", {"depth", "bases", "qualities", "alleles", "mates", "strands", "likelihoods", "hml"});

		logfile().startIndent("Will print the following pileup fields (parameter 'fields'):");
		for (const auto& field: fields) {
			auto p = Print::min;
			for (; p < Print::max; ++p) {
				if (field == printNames[p]) {
					_printSettings[p] = true;
					logfile().list(printExpl[p], " (", printNames[p], ")");
					break;
				}
			}
			if (p == Print::max) throw coretools::TUserError("Unknown field '", field, "'! Valid fields are: ", printNames);
		}
		logfile().endIndent();

		if (fields.size() < coretools::index(Print::max)) {
			logfile().startIndent("Will not print the following pileup fields (parameter 'fields'):");
			for (auto p = Print::min; p < Print::max; ++p) {
				if (!_printSettings[p]) { logfile().list(printExpl[p], " (", printNames[p], ")"); }
			}
			logfile().endIndent();
		}

		// compile header
		std::vector<std::string> header = {"chromosome", "position"};
		if (_siteTraverser.reference()) { header.push_back("reference"); }
		if (_printSettings.get<Print::Depth>()) {
			header.push_back("depth");
			if (_siteTraverser.reference()) {
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
		if (_printSettings.get<Print::Strands>()) {
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
			_printSettings.set<Print::PrintAll>(true);
			logfile().list(
			    "Will print all sites that pass filters, including those without data. (parameter 'printAll')");
		} else {
			_printSettings.set<Print::PrintAll>(false);
			logfile().list("Will print only sites with data. (use 'printAll' to print all)");
		}
	}

	constexpr coretools::TStrongArray<std::string_view, Hist> histNames{
		{"depth", "qualities", "contexts", "allelicDepth", "transitions", "prevBases"}};
	constexpr coretools::TStrongArray<std::string_view, Hist> histExpl{
		{"Sequencing depth", "Base qualities", "Base contexts", "Allelic depth", "ref to base transition",
		 "ref to base to prevBase counts"}};

	if (parameters().exists("histograms") || _onlyHistograms) {

		const auto hists = parameters().get<std::vector<std::string>>("histograms", {"depth", "qualities", "contexts"});

		logfile().startIndent("Will print the following histograms (parameter 'histograms'):");
		for (const auto& field: hists) {
			auto h = Hist::min;
			for (; h < Hist::max; ++h) {
				if (field == histNames[h]) {
					_histSettings[h] = true;
					logfile().list(histExpl[h], " (", histNames[h], ")");
					break;
				}
			}
			if (h == Hist::max) throw coretools::TUserError("Unknown field '", field, "'! Valid fields are: ", histNames);
		}
		logfile().endIndent();
		if (hists.size() < coretools::index(Hist::max)) {
			logfile().startIndent("Will not print the following pileup fields (parameter 'fields'):");
			for (auto h = Hist::min; h < Hist::max; ++h) {
				if (!_histSettings[h]) { logfile().list(histExpl[h], " (", histNames[h], ")"); }
			}
			logfile().endIndent();
		}
		if (_histSettings.get<Hist::Depth>()) {
			_outDepthHistogram.open(_siteTraverser.outputName() + "_depthPerWindow.txt.gz", {"window", "depth"});
			_outDepthPerChromosome.open(_siteTraverser.outputName() + "_depthPerChromosome.txt.gz", {"chromosome", "depth"});
		}

		if (_histSettings.get<Hist::AllelicDepth>()) {
			logfile().list("Will assemble allelic depth up to a max depth of ", _siteTraverser.upToDepth(), ". (parameter 'readUpToDepth')");
			if (_siteTraverser.upToDepth() > 100) {
				logfile().warning("Allocating count table for a max depth of ", _siteTraverser.upToDepth(), " uses a lot of memory! Use argument readUpToDepth to limit.");
			}

			_counts.resize(_siteTraverser.upToDepth());

			if (parameters().exists("includeZero")) {
				_writeEmpty = true;
				logfile().list("Will write full table, including cells with zero counts. (parameter 'includeZero')");
			} else {
				_writeEmpty = false;
				logfile().list("Will only print cells with non-zero counts. (use 'includeZero' to print all cells)");
			}
		}

		if (_histSettings.get<Hist::Transitions>()) {
			if (!_siteTraverser.reference()) {
				logfile().warning("Cannot count reference to base transitions without reference!");
				_histSettings.set<Hist::Transitions>(false);
			} else {
				using genometools::Base;
				logfile().list("Will create count table of reference-base to data-base transitions.");

				std::vector<std::string> header{"Chr", "Mate", "Strand", "End", "Pos"};
				for (auto ref = Base::min; ref < Base::max; ++ref) {
					for (auto b = Base::min; b < Base::max; ++b) {
						header.push_back(coretools::str::toString(ref, "-", b));
					}
				}
				_outTransitions.open(_siteTraverser.outputName() + "_transitions.txt.gz", header);
				_outTransitionsRel.open(_siteTraverser.outputName() + "_transitionsRel.txt.gz", header);
				_outTransitionsPsi.open(_siteTraverser.outputName() + "_transitionsPsi.txt.gz", header);
				_outTransitionsRho.open(_siteTraverser.outputName() + "_transitionsRho.txt.gz", header);
			}
		}

		if (_histSettings.get<Hist::PrevBases>()) {
			if (!_siteTraverser.reference()) {
				logfile().warning("Cannot count reference-base-previousbase without reference!");
				_histSettings.set<Hist::PrevBases>(false);
			} else {
				using genometools::Base;
				logfile().list("Will create count table of reference-base-prevousbase.");
			}
		}
	} else {
		logfile().startIndent("Will not output histograms (use 'histograms' to do so):");
		for (auto h = Hist::min; h < Hist::max; ++h) {
			logfile().list(histExpl[h], " (", histNames[h], ")");
		}
		logfile().endIndent();
	}
}

void TPileup::_addHist(const GenotypeLikelihoods::TSite& Site) {
	using genometools::Base;

	// write histograms
	if (_histSettings.get<Hist::Depth>()) {
		_depthPerSite.add(Site.depth());
		_depthPerSitePerChromosome.add(Site.depth());
	}

	if (_histSettings.get<Hist::AllelicDepth>()) {
		const auto alleleCounts = Site.countAlleles();
		_counts.addSite(alleleCounts);
	}

	if (_histSettings.get<Hist::Qualities>()) {
		for (auto &b : Site) {
			if (b.base != genometools::Base::N) { _qualDist.add(b.readGroupID, b.recalQuality.get()); }
		}
	}

	if (_histSettings.get<Hist::Contexts>()) {
		for (auto &b : Site) { _contextDist.add(b.recalQuality.get(), coretools::index(b.context())); }
	}

	if (_histSettings.get<Hist::Transitions>()) {
		for (auto &b : Site) {
			if (b.base == Base::N) continue;

			const auto p = b.dist(b.end()).pseudo();
			auto &trans  = _transitionsChr[b.mate()][b.strand()][b.end()];
			if (trans.size() <= p) { trans.resize(p + 1); }
			++trans[p][Site.refBase][b.base];
		}
	}

	if (_histSettings.get<Hist::PrevBases>()) {
		for (auto &b : Site) {
			if (b.base == Base::N) continue;
			++_prevBases[b.mate()][b.strand()][Site.refBase][b.base][b.previousBase()];
		}
	}
}

void TPileup::_printSite(const GenotypeLikelihoods::TSite &Site) {
	using genometools::Base;
	_out.write(_siteTraverser.curChr().name(),
			   _siteTraverser.position().position() + 1); // positions are zero-based internally

	if (_siteTraverser.reference()) { _out.write(Site.refBase); }
	if (_printSettings.get<Print::Depth>()) {
		_out.write(Site.depth());
		if (_siteTraverser.reference()) { _out.write(Site.refDepth(), Site.depth() - Site.refDepth()); }
	}
	if (_printSettings.get<Print::Bases>()) { _out.write(Site.getBases()); }
	if (_printSettings.get<Print::SampleBases>()) {
		const auto sBases = Site.sampleBases();
		if (sBases.empty()) {
			_out.write("-", 0, 0, 0, 0, 0);
		} else {
			coretools::TStrongArray<size_t, Base> counts{};
			size_t tot = 0;
			for (const auto b : sBases) {
				_out.writeNoDelim(b);
				++counts[b];
				++tot;
			}
			_out.writeDelim();
			_out.write(counts, tot);
		}
	}
	if (_printSettings.get<Print::Qualities>()) { _out.write(Site.getQualities()); }
	if (_printSettings.get<Print::Alleles>()) {
		const auto alleleCounts = Site.countAlleles();
		_out.write(alleleCounts[Base::A], alleleCounts[Base::C], alleleCounts[Base::G], alleleCounts[Base::T]);
		_out.write((int)coretools::numNonZero(alleleCounts));
	}
	if (_printSettings.get<Print::Mates>()) {
		const auto mateCounts = Site.countMates();
		_out.write(mateCounts);
	}
	if (_printSettings.get<Print::Strands>()) {
		const auto strandCounts = Site.countFwdRev();
		_out.write(strandCounts);
	}
	if (_printSettings.get<Print::Likelihoods>() || _printSettings.get<Print::HML>()) {
		const auto genoLik = _siteTraverser.errorModels().calculateGenotypeLikelihoods(Site);
		const auto g       = genometools::Genotype(std::max_element(genoLik.begin(), genoLik.end()) - genoLik.begin());
		if (_printSettings.get<Print::Likelihoods>()) { _out.write(genoLik, genometools::toString(g)); }
		if (_printSettings.get<Print::HML>()) { _out.write(genometools::isHeterozygous(g)); }
	}
	_out.endln();
}

void TPileup::_traverseSites() {
	logfile().startIndent("Writing pileup");

	for (; !_siteTraverser.endOfChrs(); _siteTraverser.nextChr()) {
		for (; !_siteTraverser.endOfCurChr(); _siteTraverser.nextSite()) {
			const auto &site = _siteTraverser.site();
			_addHist(site);
			if (!_onlyHistograms && (!site.empty() || _printSettings.get<Print::PrintAll>())) {
				_printSite(site);
			}
			if (_siteTraverser.winChanged()) {
				const auto& win = _siteTraverser.window();
				if (_histSettings.get<Hist::Depth>()) {
					_outDepthHistogram
						.writeNoDelim(_siteTraverser.curChr().name(), ':', win.from().position() + 1, '-',
									  win.to().position())
						.writeDelim();
					_outDepthHistogram.writeln(win.depth());
				}
			}
		}
		_endChromosome(_siteTraverser.curChr());
	}
	logfile().endIndent();
}

void TPileup::_endChromosome(const genometools::TChromosome &Chr) {
	using genometools::Base;
	using BAM::End;
	using BAM::Strand;
	using BAM::Mate;
	if (_histSettings.get<Hist::Depth>()) {
		_outDepthPerChromosome.writeln(Chr.name(), _depthPerSitePerChromosome.mean());
		_depthPerSitePerChromosome.clear();
	}
	if (_histSettings.get<Hist::Transitions>()) {
		impl::writeTransitions(_transitionsChr, Chr.name(), _outTransitions, _outTransitionsRel, _outTransitionsPsi,
						  _outTransitionsRho);
		for (auto mate = Mate::min; mate < Mate::max; ++mate) {
			for (auto strand = Strand::min; strand < Strand::max; ++strand) {
				for (auto end = End::min; end < End::max; ++end) {
					auto &chr = _transitionsChr[mate][strand][end];
					auto &tot = _transitionsTot[mate][strand][end];
					if (tot.size() < chr.size()) { tot.resize(chr.size(), {}); }
					for (auto ref = Base::min; ref < Base::max; ++ref) {
						for (auto b = Base::min; b < Base::max; ++b) {
							for (size_t i = 0; i < chr.size(); ++i) {
								tot[i][ref][b] += chr[i][ref][b];
							}
						}
					}
					chr.clear();
				}
			}
		}
	}

}

void TPileup::run() {
	using genometools::Base;
	using BAM::Mate;
	using BAM::Strand;

	_traverseSites();

	if (_histSettings.get<Hist::Depth>()) {
		// write distribution
		logfile().list("Writing depth per site distribution to file '", _siteTraverser.outputName(), "_depthPerSiteHistogram.txt.gz'");
		logfile().list("Writing average depth per window to file '", _siteTraverser.outputName(), "_depthPerWindow.txt.gz'");
		logfile().list("Writing average depth per chromosome to file '", _siteTraverser.outputName(), "_depthPerChromosome.txt.gz'");
		_depthPerSite.write(_siteTraverser.outputName() + "_depthPerSiteHistogram.txt.gz", "depth");
	}

	if (_histSettings.get<Hist::Qualities>()) {
		// print distribution
		const auto outputFileName = _siteTraverser.outputName() + "_qualHistogram.txt.gz";
		logfile().list("Writing quality distribution to '", outputFileName, "'.");
		coretools::TOutputFile out(outputFileName, {"readGroup", "quality", "counts"});

		// get read group names
		std::vector<std::string> readGroupNames;
		_siteTraverser.bamFile().readGroups().fillVectorWithNames(readGroupNames);
		// write combined
		_qualDist.resize(readGroupNames.size()); // make sure it has the right size
		_qualDist.writeCombined(out, "allReadGroups");
		_qualDist.write(out, readGroupNames);
	}

	if (_histSettings.get<Hist::Contexts>()) {
		// write counts
		const auto outputFileName = _siteTraverser.outputName() + "_contextInformation.txt.gz";
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
		const auto outputFileName = _siteTraverser.outputName() + "_allelicDepth.txt.gz";
		logfile().list("Writing allelic depth table to '", outputFileName, "' ...");
		_counts.write(outputFileName, _writeEmpty);
		logfile().done();
	}

	if (_histSettings.get<Hist::Transitions>()) {
		impl::writeTransitions(_transitionsTot, "All", _outTransitions, _outTransitionsRel, _outTransitionsPsi,
							   _outTransitionsRho);
	}

	if (_histSettings.get<Hist::PrevBases>()) {
		std::vector<std::string> header{"Mate", "Strand", "ref", "base"};
		for (auto b = Base::min; b <= Base::max; ++b) { header.push_back(toString(b)); }
		coretools::TOutputFile outPrevBase(_siteTraverser.outputName() + "_prevBases.txt.gz", header);
		coretools::TOutputFile outPrevBaseRel(_siteTraverser.outputName() + "_prevBasesRel.txt.gz", header);

		for (auto mate = Mate::min; mate < Mate::max; ++mate) {
			if (_prevBases[mate][Strand::Fwd][Base::A][Base::A][Base::A] == 0) continue; // assume mate not available
			for (auto strand = Strand::min; strand < Strand::max; ++strand) {
				for (auto ref = Base::min; ref < Base::max; ++ref) {
					for (auto b = Base::min; b < Base::max; ++b) {
						outPrevBase.write(mate, strand, ref, b);
						outPrevBaseRel.write(mate, strand, ref, b);
						size_t tot = 0;
						for (auto p = Base::min; p <= Base::max; ++p) {
							outPrevBase.write(_prevBases[mate][strand][ref][b][p]);
							tot += _prevBases[mate][strand][ref][b][p];
						}
						outPrevBase.endln();
						for (auto p = Base::min; p <= Base::max; ++p) {
							if (tot > 0) outPrevBaseRel.write(double(_prevBases[mate][strand][ref][b][p])/tot);
							else outPrevBaseRel.write(0.);
						}
						outPrevBaseRel.endln();
					}
				}
			}
		}
	}
}


} // namespace GenomeTasks
