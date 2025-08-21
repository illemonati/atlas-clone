#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include <algorithm>
#include <type_traits>

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "genometools/GenomePositions/TChromosomes.h"

#include "TBamWindows.h"
#include "TGenome.h"
#include "TWindow.h"

namespace GenomeTasks {

enum class WindowType : bool {SingleBam, MultiBam};

template<WindowType Type>
class TBamWindowTraverser {
	constexpr static bool isSingle = Type == WindowType::SingleBam;
	using GType = std::conditional_t<isSingle, TGenome, std::vector<TGenome> >;

	void _fillWindow(TGenome &genome, GenotypeLikelihoods::TWindow &Window) {
		BAM::TAlignment alignment;

		// first, use last read from last window, before reading next
		do {
			const auto curPos = genome.bamFile().curRead().position;
			const auto maxLen = genome.bamFile().curRead().cigar.lengthRead();

			if (curPos >= Window.to()) break;                                    // too far
			if (curPos.position() + maxLen < Window.from().position()) continue; // too short

			_windows.parser().fill(genome, alignment);

			if (alignment.from() < Window.to() && alignment.to() > Window.from()) {
				Window.addAlignment(alignment);
			}
		} while (genome.bamFile().readNextAlignmentThatPassesFilters());
	}

	void _fillAlignments(GenotypeLikelihoods::TWindow &Window) {
		coretools::instances::logfile().listFlushTime("Reading data ...");

		if constexpr (isSingle) {
			_fillWindow(_genome, Window);
		} else {
			for (auto& g: _genome) {
				_fillWindow(g, Window);
			}
		}
		_windows.fillSites(Window);

		coretools::instances::logfile().doneTime();

		_windows.filter(Window);
	}

	static GType _initGenome() {
		if constexpr (isSingle) {
			return {true};
		} else {
			const auto bams   = coretools::instances::parameters().get<std::vector<std::string>>("bam");
			const auto filter = BAM::TBamFilters(true);
			std::vector<TGenome> vec;
			vec.reserve(bams.size());
			vec.emplace_back(bams.front(), true, 0);
			for (size_t i = 1; i < bams.size(); ++i) {
				vec.emplace_back(bams[i], false, i);
				vec.back().bamFile().setFilters(vec.front().bamFile().filters());
			}
			return vec;
		}
	}

	static const TGenome& _front(const GType& Genome) {
		if constexpr (isSingle) {
			return Genome;
		} else {
			return Genome.front();
		}
	}

	size_t _nextRefID() {
		if constexpr (isSingle) {
			return _genome.bamFile().curPosition().refID();
		} else {
			size_t refID = -1;
			for (const auto& g : _genome) {
				refID = std::min(refID, g.bamFile().curPosition().refID());
			}
			return refID;
		}
	}

protected:
	GType _genome = _initGenome();
	TBamWindows _windows{_chromosomes(_genome)};

	static const genometools::TChromosomes& _chromosomes(const GType& Genome) {
		return _front(Genome).bamFile().chromosomes();
	}

	void _traverseBAMWindows() {
		using coretools::instances::logfile;

		if constexpr (isSingle) {
			_genome.bamFile().startProgressReporting();
			if (!_genome.bamFile().readNextAlignmentThatPassesFilters()) {
				throw coretools::TUserError("No read of file '", _genome.bamFile().filename(), "' passes filters. Are Readgroup IDs set?");
			} 

			logfile().startIndent("Traversing BAM file in windows:");
		} else {
			for (auto &g : _genome) {
				g.bamFile().startProgressReporting(false);
				if (!g.bamFile().readNextAlignmentThatPassesFilters()) {
					throw coretools::TUserError("No read of file '", g.bamFile().filename(), "' passes filters. Are Readgroup IDs set?");
				}
			}
			std::string_view file = "files";
			if (_genome.size() < 2) file.remove_suffix(1);
			logfile().startIndent("Traversing BAM ", file, " in windows:");
		}

		bool lastSkipped = false;
		for (const auto &chr : _chromosomes(_genome)) {
			_startChromosome(chr);

			if (_windows[chr.refID()].empty() || !chr.inUse() || _nextRefID() > chr.refID()) {
				if (!lastSkipped) {
					logfile().listFlush("Chromosome(s) ");
					lastSkipped = true;
				}
				logfile().flush("'", chr.name(), "',");
				_endChromosome(chr);
				continue;
			}
			if (lastSkipped) {
				logfile().write("are empty/masked.");
				lastSkipped = false;
			}

			logfile().startNumbering("Parsing chromosome '" + chr.name() + "':");

			GenotypeLikelihoods::TWindow window(chr.refID(), chr.name());
			window.allowDownsampling(_windows.allowDownsampling());
			for (const auto &gWindow : _windows[chr.refID()]) {
				window.move(gWindow, _windows.regions());
				logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ",
								 _windows[chr.refID()].size(), " on '", chr.name(), "':");

				coretools::TTimer _windowTimer;
				_fillAlignments(window);
				if (window.passedFilters()) {
					_handleWindow(window);
					logfile().list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
				}
			}
			_endChromosome(chr);
			logfile().endNumbering();
		}
		if (lastSkipped) {
			logfile().write("are empty/masked.");
		}
		logfile().endIndent();

		if constexpr (isSingle) {
			_genome.bamFile().printEndWithSummary(_genome.outputName());
		} else {
			for (auto& g: _genome) {
				g.bamFile().printEndWithSummary(g.outputName(), false);
			}
		}
	}

	virtual void _handleWindow(GenotypeLikelihoods::TWindow &window)   = 0;
	virtual void _startChromosome(const genometools::TChromosome &Chr) = 0;
	virtual void _endChromosome(const genometools::TChromosome &Chr)   = 0;

public:
	virtual ~TBamWindowTraverser() = default;
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
