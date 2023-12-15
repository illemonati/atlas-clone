#include "TAlignmentList.h"
#include "coretools/Files/TFile.h"
namespace BAM{

TAlignmentList::TAlignmentList(std::string_view filename){
	addFromFile(filename);
};

void TAlignmentList::addFromFile(std::string_view filename){
	coretools::TInputFile in(filename, 1);
	std::vector<std::string> vec;
	while(in.read(vec)){
		add(vec[0]);
	}
};

void TAlignmentList::add(std::string_view alignment){
	_list.emplace(alignment);
};

void TAlignmentList::remove(const std::string& alignment){
	_list.erase(alignment);
};

void TAlignmentList::clear(){
	_list.clear();
};

bool TAlignmentList::isInBlacklist(const std::string & alignment) const{
	return _list.find(alignment) != _list.end();
};
}
