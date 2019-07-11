#ifndef BAMREADERSEQLIB_H
#define BAMREADERSEQLIB_H

#include "IOAbstractClasses/SamHeader.h"
#include "SamHeaderSeqLib.h"
#include "IOAbstractClasses/BamAlignment.h"
#include "BamAlignmentSeqLib.h"
#include "IOAbstractClasses/RefVector.h"
#include "RefVectorSeqLib.h"
#include "IOAbstractClasses/BamReader.h"

#include "SeqLib/BamReader.h"
#include "SeqLib/GenomicRegion.h"

class BamReaderSeqLib : public BamReader
{
private:
    SeqLib::BamReader bamReader;
    RefVectorSeqLib refVectors;
    BamAlignmentSeqLib nextAlignment;
    SamHeaderSeqLib samHeader;
public:
    BamReaderSeqLib();

    bool Open(std::string& filepath);
    bool Close();
    bool LocateIndex();
    SamHeader* GetHeader();
    bool Jump(int refid, int pos);
    bool GetNextAlignment(BamAlignment*& alignment);
    int GetReferenceCount();
    RefVector& GetReferenceData();
    bool CreateIndex();
    void Rewind();
    int tell();

    ~BamReaderSeqLib(){}
};

#endif // BAMREADERSEQLIB_H
