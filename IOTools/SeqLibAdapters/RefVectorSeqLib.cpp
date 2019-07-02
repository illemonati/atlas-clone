#include "RefVectorSeqLib.h"

RefVectorSeqLib::RefVectorSeqLib(){}

RefVectorSeqLib::RefVectorSeqLib(SeqLib::HeaderSequenceVector references)
{
    this->references=references;
}

void RefVectorSeqLib::push_back(RefData & refdata)
{
    RefDataSeqLib child = dynamic_cast<RefDataSeqLib&>(refdata);
    references.push_back(child.reference);
}

void RefVectorSeqLib::clear()
{
    references.clear();
}

SeqLib::HeaderSequenceVector RefVectorSeqLib::GetRefVector()
{
    return references;
}


RefDataSeqLib::RefDataSeqLib(std::string & name, uint32_t length)
{
    reference = SeqLib::HeaderSequence(name,length);
}
