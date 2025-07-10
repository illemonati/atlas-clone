#include "TAlignmentList.h"

#include "coretools/Files/TLineReader.h"

namespace BAM{

void TAlignmentList::addFromFile(std::string_view filename){
	std::vector<std::string> vec;
	for (coretools::TLineReader in(filename); !in.empty(); in.popFront()) {
		add(in.front());
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
