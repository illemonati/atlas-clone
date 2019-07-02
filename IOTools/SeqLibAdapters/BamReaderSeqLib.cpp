#include "BamReaderSeqLib.h"

BamReaderSeqLib::BamReaderSeqLib(){}

bool BamReaderSeqLib::Open(std::string &filepath)
{
    return bamReader.Open(filepath);;
}

bool BamReaderSeqLib::Close()
{
    return bamReader.Close();
}

bool BamReaderSeqLib::LocateIndex()
{
    return true;
}

SamHeader* BamReaderSeqLib::GetHeader()
{
    samHeader = SamHeaderSeqLib(bamReader.Header());  
    return &samHeader;
}

bool BamReaderSeqLib::Jump(int refid, int pos)
{
    SeqLib::GenomicRegion region = SeqLib::GenomicRegion(refid,pos,bamReader.Header().GetSequenceLength(refid));
    return bamReader.SetRegion(region);
}

bool BamReaderSeqLib::GetNextAlignment(BamAlignment*& alignment)
{
    nb_call++;
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    SeqLib::BamRecord record;
    bool success = bamReader.GetNextRecord(record);
    nextAlignment = BamAlignmentSeqLib(record);
    alignment = &nextAlignment;
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    duration += std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    return success;
}

int BamReaderSeqLib::GetReferenceCount()
{
    return bamReader.Header().NumSequences();
}

RefVector &BamReaderSeqLib::GetReferenceData()
{
    refVectors = RefVectorSeqLib(bamReader.Header().GetHeaderSequenceVector());
    return refVectors;
}

bool BamReaderSeqLib::CreateIndex()
{
    //SeqLib auto create the index on openning if it doesn't exist
    return true;
}

void BamReaderSeqLib::Rewind()
{
    bamReader.Reset();
}

int BamReaderSeqLib::tell()
{
    return bamReader.Tell();
}
