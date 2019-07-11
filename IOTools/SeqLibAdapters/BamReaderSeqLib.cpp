#include "BamReaderSeqLib.h"

BamReaderSeqLib::BamReaderSeqLib(){}

bool BamReaderSeqLib::Open(std::string &filepath)
{
    bool res = bamReader.Open(filepath);
    std::cout << "RESULT=" << bamReader.Header().AsString() << "\n";
    return res;
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
    static int count=0;
    std::cout << "Number alignment="<<std::to_string(count)<<"\n";
    count++;
    SeqLib::BamRecord record;
    bool success = bamReader.GetNextRecord(record);
    nextAlignment = BamAlignmentSeqLib(record);
    alignment = &nextAlignment;
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
