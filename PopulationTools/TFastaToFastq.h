#ifndef TFASTATOFASTQ_H_
#define TFASTATOFASTQ_H_

#include "TBamTraverser.h"


namespace PopulationTools {

struct TFastaToFastq {
	void run();
}

class TReadCenterAlignedBase : public TBamReadTraverser<ReadType::Parsed> {
	private: 
		void _handleAlignment(BAM::TAlignment& alignment);

	public:
	TReadCenterAlignedBase();
	void run();
};
} // namespace PopulationTools

#endif
