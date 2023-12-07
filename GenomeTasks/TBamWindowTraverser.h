#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include "TBamWindows.h"
#include "TGenome.h"
#include "TWindow.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/splitters.h"
#include "genometools/GenomePositions/TChromosomes.h"
#include <type_traits>

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
			const auto curPos = genome.bamFile().curPosition();
			const auto maxLen = genome.bamFile().curCIGAR().lengthSequenced();

			if (curPos >= Window.to()) break;                                    // too far
			if (curPos.position() + maxLen < Window.from().position()) continue; // too short

			_windows.parser().fill(genome, alignment);

			if (alignment.lastAlignedPositionWithRespectToRef() >= Window.from()) { Window.addAlignment(alignment); }
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
			return {BAM::TBamFilters(true)};
		} else {
			const auto bams   = coretools::instances::parameters().get<std::vector<std::string>>("bam");
			const auto filter = BAM::TBamFilters(true);
			std::vector<TGenome> vec;
			vec.reserve(bams.size());
			if (bams.size() == 1) {
				vec.emplace_back(bams.front(), filter);
			} else {
				for (size_t i = 0; i < bams.size(); ++i) {
					vec.emplace_back(bams[i], filter, i);
				}
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

	static const genometools::TChromosomes& _chromosomes(const GType& Genome) {
		return _front(Genome).bamFile().chromosomes();
	}

protected:
	GType _genome = _initGenome();
	TBamWindows _windows{_chromosomes(_genome)};

	void _traverseBAMWindows() {
		using coretools::instances::logfile;

		if constexpr (isSingle) {
			_genome.bamFile().startProgressReporting();
			_genome.bamFile().readNextAlignmentThatPassesFilters();
			logfile().startIndent("Traversing BAM file in windows:");
		} else {
			for (auto &g : _genome) {
				g.bamFile().startProgressReporting(false);
				g.bamFile().readNextAlignmentThatPassesFilters();
			}
			std::string_view file = "files";
			if (_genome.size() < 2) file.remove_suffix(1);
			logfile().startIndent("Traversing BAM ", file, " in windows:");
		}

		for (const auto &chr : _chromosomes(_genome)) {
			if (!chr.inUse()) continue;

			logfile().startNumbering("Parsing chromosome '" + chr.name() + "':");
			_startChromosome(chr);

			GenotypeLikelihoods::TWindow window(chr.name());
			for (const auto &gWindow : _windows[chr.refID()]) {
				window.move(gWindow);
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
	TBamWindowTraverser() {}
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
