#include "TAlignmentList.h"
#include "coretools/Files/TInputFile.h"

namespace BAM{

void TAlignmentList::addFromFile(std::string_view filename){
	coretools::TInputFile in(filename, coretools::FileType::NoHeader);
	std::vector<std::string> vec;
	for (coretools::TInputFile in(filename, coretools::FileType::NoHeader); !in.empty(); in.popFront()) {
		add(in.get(0));
	}
}

void TAlignmentList::add(std::string_view alignment){
	_list.emplace(alignment);
}

void TAlignmentList::remove(const std::string& alignment){
	_list.erase(alignment);
}

void TAlignmentList::clear(){
	_list.clear();
}

bool TAlignmentList::isInBlacklist(const std::string & alignment) const{
	return _list.find(alignment) != _list.end();
}

} // namespace BAM
