#ifndef GENOMETASKS_TFROMTO_H_
#define GENOMETASKS_TFROMTO_H_

#include "TSiteTraverser.h"
#include "coretools/Files/TOutputFile.h"

namespace GenomeTasks{

class TFromTo {
private:
	BAM::TSiteTraverser _siteTraverser;

	coretools::TOutputFile _out;
	void _traverseSites();
public:
	TFromTo();
	void run();
};

}; // end namespace



#endif /* GENOMETASKS_TDEPTHWRITER_H_ */
