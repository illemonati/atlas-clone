#ifndef TBAMWINDOWTRAVERSER_H_
#define TBAMWINDOWTRAVERSER_H_

#include "TBamWindows.h"
#include "TGenome.h"
#include "TWindow.h"
#include "coretools/Main/TLog.h"
#include <type_traits>

namespace GenomeTasks {

enum class WindowType : bool {SingleBam, MultiBam};

template<WindowType Type>
class TBamWindowTraverser {
	using GType = std::conditional_t<Type == WindowType::SingleBam, TGenome, std::vector<TGenome> >;

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

		_fillWindow(_genome, Window);
		_windows.fillSites(Window);

		coretools::instances::logfile().doneTime();

		_windows.filter(Window);
	}

	static GType _initGenome() {
		return GType{BAM::TBamFilters(true)};
	}

protected:
	GType _genome;
	TBamWindows _windows;

	void _traverseBAMWindows() {
		using coretools::instances::logfile;

		_genome.bamFile().startProgressReporting();
		_genome.bamFile().readNextAlignmentThatPassesFilters();
		logfile().startIndent("Traversing BAM file in windows:");

		for (const auto &chr : _genome.bamFile().chromosomes()) {
			if (!chr.inUse()) continue;

			logfile().startNumbering("Parsing chromosome '" + chr.name() + "':");
			_handleChromosome(chr);

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
			logfile().endIndent();
		}
		logfile().endIndent();

		_genome.bamFile().printEndWithSummary(_genome.outputName());
	}

	virtual void _handleWindow(GenotypeLikelihoods::TWindow &window)    = 0;
	virtual void _handleChromosome(const genometools::TChromosome &Chr) = 0;

public:
	TBamWindowTraverser() : _genome(_initGenome()), _windows(_genome.bamFile().chromosomes()) {}
};

} // namespace GenomeTasks

#endif /* GENOME_H_ */
