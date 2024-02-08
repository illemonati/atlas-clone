#include "TBaseFilter.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/concatenateString.h"
namespace GenomeTasks {

using coretools::instances::parameters;
using coretools::instances::logfile;
//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
TQualityFilter::TQualityFilter() {
	if(parameters().exists("filterBaseQual")){
		parameters().fill("filterBaseQual", _range);
		if (_range.within(genometools::PhredIntProbability(0))){ UERROR("Base quality filter of 0 is not allowed (parameter 'filterBaseQual')"); }
		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
	} else {
		_range.set(genometools::PhredIntProbability(1), true, genometools::PhredIntProbability(93), true);
		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + ". (use 'filterBaseQual' to change)");
	}
};

//-------------------------------------
// TContextFilter
//-------------------------------------
TContextFilter::TContextFilter(){
	using namespace genometools;

	_keptContexts.fill(true);
	_filter = false;
	if(parameters().exists("ignoreContexts")){
		std::vector<std::string> contexts;
		parameters().fill("ignoreContexts", contexts);

		if(contexts.size() > 0){
			for(auto& c : contexts){
				if(c.size() != 2){
					UERROR("Context ", c, " does not consist of two bases! (parameter 'ignoreContexts')");
				}

				const Base first  = char2base(c[0]);
				const Base second = char2base(c[1]);

				if(base2char(first) != c[0] || base2char(second) != c[1]){
					UERROR("Unable to understand context '", c, "'!  (parameter 'ignoreContexts')");
				}

				//save context
				_keptContexts[baseContext(first, second)] = false;
			}

			std::vector<std::string> rep;
			for(auto i = BaseContext::min; i <= BaseContext::max; ++i){ //including max
				if(!_keptContexts[i]){
					rep.push_back(toString(i));
				}
			}
			logfile().list("Will ignore the following contexts: " + coretools::str::concatenateString(rep, ", ")  + ". (parameter 'ignoreContexts')");
			_filter = true;
		}
	}

	if(!_filter){
		logfile().list("Will keep bases regardless of base context. (use 'ignoreContexts' to filter)");
	}
};
}
