#include "TBamWindowTraverser.h"

#include <algorithm>
#include <filesystem>

#include "TWindow.h"
#include "coretools/Main/TLog.h"

namespace GenomeTasks {

using coretools::instances::logfile;

void TBamWindowTraverser::_fillAlignments(GenotypeLikelihoods::TWindow &Window) {
	// measure runtime
	logfile().listFlushTime("Reading data ...");
	BAM::TAlignment alignment;

	// first, use last read from last window, before reading next
	do {
		const auto curPos = _genome.bamFile().curPosition();
		const auto maxLen = _genome.bamFile().curCIGAR().lengthSequenced();
		
		if (curPos >= Window.to()) break; // too far
		if (curPos.position() + maxLen < Window.from().position()) continue; // too short

		_windows.parser().fill(_genome, alignment);

		if (alignment.lastAlignedPositionWithRespectToRef() >= Window.from()) {
			Window.addAlignment(alignment);
		}
	} while (_genome.bamFile().readNextAlignmentThatPassesFilters());

	// fill sites
	_windows.fillSites(Window);

	logfile().doneTime();

	// apply filters
	_windows.filter(Window);
}

void TBamWindowTraverser::_traverseBAMWindows() {
	using GenotypeLikelihoods::TWindow;

	_genome.bamFile().startProgressReporting();
	logfile().startIndent("Traversing BAM file in windows:");
	_genome.bamFile().readNextAlignmentThatPassesFilters();
	for (const auto &chr : _genome.bamFile().chromosomes()) {
		if (!chr.inUse()) continue;
		logfile().startNumbering("Parsing chromosome '" + chr.name() + "':");
		_handleChromosome(chr);
		TWindow window(chr.name());

		for (const auto& gWindow: _windows[chr.refID()]) {
			window.move(gWindow);
			logfile().number("Window [", window.from().position() + 1, ", ", window.to().position(), "] of ",
							 _windows[chr.refID()].size(), " on '", chr.name(), "':");

			coretools::TTimer _windowTimer;
			_fillAlignments(window);
			if (window.passedFilters()) {
				// do stuff in derived classes
				_handleWindow(window);

				// report end of window calculations
				logfile().list("Total computation time for this window was ", _windowTimer.formattedTime(), ".");
			}
		}
		logfile().endIndent();
	}
	logfile().endIndent();

	// report
	_genome.bamFile().printEndWithSummary(_genome.outputName());
}

} // namespace GenomeTasks
