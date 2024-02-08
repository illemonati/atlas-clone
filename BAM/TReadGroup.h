
#ifndef TREADGROUP_H_
#define TREADGROUP_H_

#include <cstddef>
#include <string>

namespace BAM {

//---------------------------------------------------------------
// TReadGroup
//---------------------------------------------------------------
struct TReadGroup {
public:
	static constexpr size_t noID = -1;

	size_t id           = -1; // internal ID
	std::string name_ID = "No Read Group";
	std::string barcodeSequence_BC;
	std::string sequencingCenter_CN;
	std::string description_DS;
	std::string productionDate_DT;
	std::string flowOrder_FO;
	std::string keySequence_KS;
	std::string library_LB;
	std::string program_PG;
	std::string predictedInsertSize_PI;
	std::string sequencingTechnology_PL;
	std::string platformModel_PM;
	std::string platformUnit_PU;
	std::string sample_SM;

	// flags
	bool inUse         = false; // read groups not in use are ignored when reading
	bool writeToHeader = false; // is false if read group is not in use or replaced by new one

	TReadGroup() = default;
	TReadGroup(size_t ID, std::string_view Name) : id(ID), name_ID(Name), inUse(true), writeToHeader(true){};

	std::string compileSamHeader() const;

	bool operator<(const TReadGroup &rhs) const { return name_ID < rhs.name_ID; };
	bool operator<(std::string_view rhs) const { return name_ID < rhs; };
	bool operator==(std::string_view name) const { return name_ID == name; };
	friend bool operator<(std::string_view lhs, const TReadGroup &rhs) { return lhs < rhs.name_ID; };
};

} // namespace BAM

#endif
