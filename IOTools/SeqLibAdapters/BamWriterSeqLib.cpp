#include "BamWriterSeqLib.h"


bool BamWriterSeqLib::Open(std::string &filepath, SamHeader &header, RefVector &references)
{
    bool success = bamWriter.Open(filepath);

    RefVectorSeqLib sequences = dynamic_cast<RefVectorSeqLib&>(references);
    auto bamHeader=SeqLib::BamHeader(sequences.GetRefVector());

    SamHeaderSeqLib samHeader= dynamic_cast<SamHeaderSeqLib&>(header);
    samHeader.SetSamHeader(bamHeader);

    bamWriter.SetHeader(bamHeader);
    bamWriter.WriteHeader();
    return success;
}

bool BamWriterSeqLib::Close()
{
    bool res= bamWriter.Close();
    bamWriter.BuildIndex();
    return res;
}

bool BamWriterSeqLib::SaveAlignment(BamAlignment &alignment)
{
    BamAlignmentSeqLib & child= dynamic_cast<BamAlignmentSeqLib&>(alignment);
    return bamWriter.WriteRecord(child.GetSeqLibRecord());
}
