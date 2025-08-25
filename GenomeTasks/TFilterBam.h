#ifndef TFILTERBAM_H_
#define TFILTERBAM_H_

#include "TSequencedData.h"
#include "TWaitingListBamTraverser.h"

namespace BAM {
class TAlignment;
}
namespace GenomeTasks {

//-----------------------------------------
// TBamFilter
//-----------------------------------------
class TFilterBam final:public BAM::TWaitingListBamTraverser {
private:
	size_t _makeSingle = 0;

	void _handleMates(BAM::TWaitingAlignment &lhs, BAM::TWaitingAlignment &rhs) override;
	void _handleSingle(BAM::TWaitingAlignment &lhs) override;
	void _handleOrphan(BAM::TWaitingAlignment &lhs) override;

public:
	TFilterBam();
	void run() {
		traverseBAM();
	}
};

} // namespace GenomeTasks::BamFilter

#endif /* TBAMFILTER_H_ */
