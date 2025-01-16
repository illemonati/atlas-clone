/*
 * TAlignmentMerger.h
 *
 *  Created on: Oct 20, 2022
 *      Author: phaentu
 */

#ifndef GENOMETASKS_TALIGNMENTMERGER_H_
#define GENOMETASKS_TALIGNMENTMERGER_H_


#include <cstddef>
#include <utility>
namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }

namespace GenomeTasks::AlignmentMerger {


//-----------------------------------------
// TAlignmentMergerReadGroupSettings
//-----------------------------------------
class TBase{
private:
	size_t overlapLengthAndMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
protected:
	size_t _merge(BAM::TAlignment & alignment, BAM::TAlignment & mate);
public:
	virtual ~TBase() = default;
	virtual size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) = 0;
	static size_t determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate);
};

class TRandomRead final : public TBase {
public:
	TRandomRead();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

class TMiddle final :public TBase{
public:
	TMiddle();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
	std::pair<size_t,bool> determineOverlapLength(const BAM::TAlignment & alignment, const BAM::TAlignment & mate);
	void sameDirectionMerge(BAM::TAlignment & alignment, BAM::TAlignment & mate, std::pair<size_t,bool> overlapLength);
};

class TFirstMate final : public TBase {
public:
	TFirstMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

class TSecondMate final : public TBase {
public:
	TSecondMate();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

class THighestQuality final : public TBase {
private:
	size_t _overlapLengthAndMergeBetterQuality(BAM::TAlignment & alignment, BAM::TAlignment & mate);
public:
	THighestQuality();
	size_t merge(BAM::TAlignment & alignment, BAM::TAlignment & mate) override;
};

} // namespace GenomeTasks::AlignmentMerger

#endif /* GENOMETASKS_TALIGNMENTMERGER_H_ */
