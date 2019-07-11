#include "BamReaderBamTools.h"

BamReaderBamTools::BamReaderBamTools(){}

bool BamReaderBamTools::Open(std::string &filepath){
    return bamReader.Open(filepath);
}

bool BamReaderBamTools::Close()
{
    return bamReader.Close();
}

bool BamReaderBamTools::LocateIndex(){
    return bamReader.LocateIndex();
}

SamHeader* BamReaderBamTools::GetHeader(){
    samHeader = SamHeaderBamTools(bamReader.GetHeader());
    return &samHeader;
}

bool BamReaderBamTools::Jump(int refid, int pos){
    return bamReader.Jump(refid,pos);
}

bool BamReaderBamTools::GetNextAlignment(BamAlignment *& alignment){
    BamTools::BamAlignment record;
    bool success = bamReader.GetNextAlignment(record);
    nextAlignment = BamAlignmentBamTools(record);
    alignment=&nextAlignment;
    return success;
}

int BamReaderBamTools::GetReferenceCount(){
    return bamReader.GetReferenceCount();
}

RefVector &BamReaderBamTools::GetReferenceData(){
    refVectors=RefVectorBamTools(bamReader.GetReferenceData());
    return refVectors;
}

bool BamReaderBamTools::CreateIndex()
{
    return bamReader.CreateIndex(BamTools::BamIndex::STANDARD);
}

void BamReaderBamTools::Rewind(){
    bamReader.Rewind();
}

int BamReaderBamTools::tell(){
    return bamReader.tell();
}
