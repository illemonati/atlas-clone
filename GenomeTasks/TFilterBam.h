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
class TFilterBam final:public TWaitingListBamTraverser {
private:
	size_t _makeSingle = 0;

	void _handleMates(TWaitingAlignment &lhs, TWaitingAlignment &rhs) override;
	void _handleSingle(TWaitingAlignment &lhs) override;
	void _handleOrphan(TWaitingAlignment &lhs) override;

public:
	TFilterBam();
	void run() {
		traverseBAM();
	}
};

} // namespace GenomeTasks::BamFilter

#endif /* TBAMFILTER_H_ */
