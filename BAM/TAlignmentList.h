
#ifndef BAM_TALIGNMENTLIST_H
#define BAM_TALIGNMENTLIST_H

#include <set>
#include <string>

namespace BAM{

class TAlignmentList{
private:
	std::set<std::string> _list;

public:
	TAlignmentList() = default;
	TAlignmentList(std::string_view filename);
	void addFromFile(std::string_view filename);
	void add(std::string_view alignment);
	void remove(const std::string& alignment);
	void clear();
	bool isInBlacklist(const std::string& alignment) const;
};
}
#endif
