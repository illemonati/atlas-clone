#ifndef REFVECTORSEQLIB_H
#define REFVECTORSEQLIB_H

#include "SeqLib/BamHeader.h"
#include "IOAbstractClasses/RefVector.h"


class RefVectorSeqLib: public RefVector
{
private:
    SeqLib::HeaderSequenceVector references;
public:
    RefVectorSeqLib();
    RefVectorSeqLib(SeqLib::HeaderSequenceVector);

    void push_back(RefData&);
    void clear();

    SeqLib::HeaderSequenceVector GetRefVector();

    ~RefVectorSeqLib(){}
};

class RefDataSeqLib : public RefData
{
private:
    SeqLib::HeaderSequence reference = SeqLib::HeaderSequence("",0);
public:
    RefDataSeqLib(std::string&, uint32_t);

    ~RefDataSeqLib(){}

    friend class RefVectorSeqLib;
};

#endif // REFVECTORSEQLIB_H
