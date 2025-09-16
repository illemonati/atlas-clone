#ifndef GENOMETASKS_TWINDOWSTATS_H_
#define GENOMETASKS_TWINDOWSTATS_H_

#include "TSiteTraverser.h"

namespace GenomeTasks {
class TWindowStats {
	BAM::TSiteTraverser _siteTraverser;
	coretools::TOutputFile _out;
	size_t _genomeWide = false;

	void _traverseWindows();
	void _traverseGenome();
	
public:
	TWindowStats();
	void run();
};
} // namespace GenomeTasks
#endif
