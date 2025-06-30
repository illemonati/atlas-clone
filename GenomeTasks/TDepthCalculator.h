#ifndef TDEPTHCALCULATOR_H_
#define TDEPTHCALCULATOR_H_

#include "TGenome.h"
#include "coretools/Containers/TView.h"
#include "coretools/Main/TParameters.h"
#include <string>

namespace GenomeTasks {

class TDepthCalculator {
	GenomeTasks::TGenome _genome{true};
public:
	void run();
};
}
#endif
