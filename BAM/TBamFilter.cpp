/*
 * TBamFilter.cpp
 *
 *  Created on: May 24, 2020
 *      Author: phaentu
 */


#include "TBamFilter.h"
#include "GenotypeTypes.h"

namespace BAM{

//-----------------------------------------------------
//TBamFileLog
//-----------------------------------------------------
TBamFileLog::TBamFileLog(const std::string filename){
	_log.open(filename);
	//_log.writeHeader({"Alignment", "isSecondMate", "Reason"});
	_log.noHeader(3); //write no header so it can be used as blacklist
};

void TBamFileLog::write(const std::string & alignmentName, const bool & isReverseStrand, const std::string & reason){
	_log << alignmentName << isReverseStrand << reason << std::endl;
};

//-----------------------------------------------------
//TBamFileFilter_base
//-----------------------------------------------------
TBamFileFilter::TBamFileFilter(){
	_counter = 0;
	_keep = true;
	_updateLog = false;
	_log = nullptr;
};

void TBamFileFilter::filterOut(const std::string & alignmentName, const bool & isReverseStrand){
	++_counter;
	if(_updateLog){
		_log->write(alignmentName, isReverseStrand, _reason);
	}
};

void TBamFileFilter::keep(){
	_keep = true;
};

void TBamFileFilter::setReason(const std::string reason){
	_reason = reason;
};

void TBamFileFilter::setLog(TBamFileLog* Log){
	_log = Log;
	_updateLog = true;
};

void TBamFileFilter::summary(TLog* logfile, uint64_t total){
	if(!_keep && _counter  > 0){
		logfile->list(_reason + ": ", _counter, " (" + coretools::str::toPercentString(_counter, total, 3) + "%)");
	}
};

//-----------------------------------------------------
//TBamFileFilterBool
//-----------------------------------------------------
void TBamFileFilterBool::filter(const std::string Reason){
	_keep = false;
	_reason = Reason;
};

bool TBamFileFilterBool::pass(const bool state, const std::string & alignmentName, const bool & isReverseStrand){
	if(state && !_keep){
		filterOut(alignmentName, isReverseStrand);
		return false;
	}
	return true;
};


//---------------------------------------------------------------
//TQualityFilter
//---------------------------------------------------------------
void TQualityFilter::_default(){
	//default values according to SAM specifications
	_filter = false;
	_range.set(genometools::PhredIntProbability(1), true, genometools::PhredIntProbability(93), true);
};

void TQualityFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	if(params.parameterExists("filterBaseQual")){
		params.fillParameter("filterBaseQual", _range);
		logfile->list("Will filter out bases with quality outside the range " + _range.rangeString() + " (parameter 'filterBaseQual')");
		_filter = true;
	} else {
		_default();
		logfile->list("Will keep bases regardless of quality (use 'filterBaseQual' to filter)");
	}
};

//-------------------------------------
// TContextFilter
//-------------------------------------
void TContextFilter::set(coretools::TParameters & params, coretools::TLog* logfile){
	using namespace genometools;
	_filter = false;
	if(params.parameterExists("ignoreContexts")){
		std::vector<std::string> contexts;
		params.fillParameterIntoContainer("ignoreContexts", contexts, ',');

		if(contexts.size() > 0){
			for(auto& c : contexts){
				if(c.size() != 2){
					throw "Context " + c + " does not consist of two bases! (parameter 'ignoreContexts')";
				}

				const Base first  = char2base(c[0]);
				const Base second = char2base(c[1]);

				if((char) first != c[0] || (char) second != c[1]){
					throw "Unable to understand context '" + c + "'!  (parameter 'ignoreContexts')";
				}

				//save context
				_keptContexts[index(baseContext(first, second))] = false;
			}

			std::vector<std::string> rep;
			for(auto i = 0; i <= index(BaseContext::NN); ++i){
				if(!_keptContexts[i]){
					rep.push_back(toString(BaseContext(i)));
				}
			}
			logfile->list("Will ignore the following contexts: " + coretools::str::concatenateString(rep, ", ")  + ". (parameter 'ignoreContexts')");
			_filter = true;
		}
	}

	if(!_filter){
		logfile->list("Will keep bases regardless of base context. (use 'ignoreContexts' to filter)");
	}
};

bool TContextFilter::pass(const TSequencedBase & base) const{
	return _keptContexts[index(base.context)];
};

//-----------------------------------------------------
//TAlignmentBlacklist
//-----------------------------------------------------
TAlignmentList::TAlignmentList(const std::string filename){
	addFromFile(filename);
};

void TAlignmentList::addFromFile(const std::string filename){
	coretools::TInputFile in(filename, 1);
	std::vector<std::string> vec;
	while(in.read(vec)){
		add(vec[0]);
	}
};

void TAlignmentList::add(const std::string & alignment){
	_list.insert(alignment);
};

void TAlignmentList::remove(const std::string & alignment){
	_list.erase(alignment);
};

void TAlignmentList::clear(){
	_list.clear();
};

bool TAlignmentList::isInBlacklist(const std::string & alignment) const{
	return _list.find(alignment) != _list.end();
};

}; //end namespace
