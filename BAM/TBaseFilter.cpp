#include "TBaseFilter.h"

#include <sys/types.h>

#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Strings/concatenateString.h"

#include "genometools/Genotypes/Base.h"
#include "genometools/Genotypes/BaseContext.h"


namespace BAM {

using coretools::instances::parameters;
using coretools::instances::logfile;
using coretools::user_assert;

//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
TQualityFilter::TQualityFilter() {
	if(parameters().exists("filterBaseQual")){
		parameters().fill("filterBaseQual", _range);
		user_assert(!_range.within(coretools::PhredInt::highest()), "Base quality filter of 0 is not allowed (parameter 'filterBaseQual')");

		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
	} else {
		_range.set(coretools::PhredInt(1), true, coretools::PhredInt(93), true);
		logfile().list("Will filter out bases with quality outside the range " + _range.rangeString() + ". (use 'filterBaseQual' to change)");
	}
};

//-------------------------------------
// TContextFilter
//-------------------------------------
TContextFilter::TContextFilter(){
	using genometools::Base;
	using genometools::BaseContext;
	using genometools::char2base;

	_keptContexts.fill(true);
	_filter = false;
	if(parameters().exists("ignoreContexts")){
		std::vector<std::string> contexts;
		parameters().fill("ignoreContexts", contexts);

		if(contexts.size() > 0){
			for(auto& c : contexts){
				user_assert(c.size() == 2, "Context ", c, " does not consist of two bases! (parameter 'ignoreContexts')");

				const Base first  = char2base(c[0]);
				const Base second = char2base(c[1]);

				user_assert(base2char(first) == c[0] && base2char(second) == c[1], "Unable to understand context '", c,
							"'!  (parameter 'ignoreContexts')");

				// save context
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
