#ifndef BAMWRITERSEQLIB_H
#define BAMWRITERSEQLIB_H

#include <string>
#include "../IOAbstractClasses/BamWriter.h"
#include "../IOAbstractClasses/SamHeader.h"
#include "SamHeaderSeqLib.h"
#include "../IOAbstractClasses/RefVector.h"
#include "RefVectorSeqLib.h"
#include "../IOAbstractClasses/BamAlignment.h"
#include "BamAlignmentSeqLib.h"

#include "SeqLib/BamWriter.h"
#include "SeqLib/BamHeader.h"

class BamWriterSeqLib : public BamWriter
{
private:
    SeqLib::BamWriter bamWriter;
public:
    bool Open(std::string& filepath, SamHeader& header, RefVector& references);
    bool Close();
    bool SaveAlignment(BamAlignment& alignment);

    ~BamWriterSeqLib(){}
};

#endif // BAMWRITERSEQLIB_H
