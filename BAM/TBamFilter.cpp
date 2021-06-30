/*
 * TBamFilter.cpp
 *
 *  Created on: May 24, 2020
 *      Author: phaentu
 */


#include "TBamFilter.h"

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

void TBamFileFilter::summary(TLog* logfile, const uint64_t & total){
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
