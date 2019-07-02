#include "BamAlignmentSeqLib.h"

BamAlignmentSeqLib::BamAlignmentSeqLib(){}

BamAlignmentSeqLib::BamAlignmentSeqLib(SeqLib::BamRecord alignment)
{
    this->alignment = alignment;
}

int BamAlignmentSeqLib::GetRefID()
{
    return alignment.ChrID();
}

int BamAlignmentSeqLib::GetMapQuality()
{
    return alignment.MapQuality();
}

int BamAlignmentSeqLib::GetLength()
{
    return alignment.Length();
}

int BamAlignmentSeqLib::GetPosition()
{
    return alignment.Position();
}

int BamAlignmentSeqLib::GetMateRefID()
{
    return alignment.MateChrID();
}

int BamAlignmentSeqLib::GetMatePosition()
{
    return alignment.MatePosition();
}

int BamAlignmentSeqLib::GetInsertSize()
{
    return alignment.InsertSize();
}

std::string BamAlignmentSeqLib::GetAlignedBases()
{
    return alignment.AlignedBases();
}

std::string BamAlignmentSeqLib::GetQueryBases()
{
    if(queryBases.empty()){
        queryBases=alignment.Sequence();
    }

    return queryBases;
}

std::string BamAlignmentSeqLib::GetQualities()
{
    if(qualities.empty()){
        qualities=alignment.Qualities();
    }
    return qualities;
}

std::string BamAlignmentSeqLib::GetName()
{
    return alignment.Qname();
}

bool BamAlignmentSeqLib::IsReverseStrand()
{
    return alignment.ReverseFlag();
}

bool BamAlignmentSeqLib::IsPaired()
{
    return alignment.PairedFlag();
}

bool BamAlignmentSeqLib::IsSecondMate()
{
    return alignment.SecondaryFlag();
}

bool BamAlignmentSeqLib::IsProperPair()
{
    return alignment.ProperPair();
}

bool BamAlignmentSeqLib::IsMapped()
{
    return alignment.MappedFlag();
}

bool BamAlignmentSeqLib::IsFailedQC()
{
    return alignment.QCFailFlag();
}

bool BamAlignmentSeqLib::IsPrimaryAlignment()
{
    return !alignment.SecondaryFlag();
}

bool BamAlignmentSeqLib::IsSupplementary()
{
    return alignment.SupplementaryFlag();
}

bool BamAlignmentSeqLib::IsDuplicate()
{
    return alignment.DuplicateFlag();
}

Cigar& BamAlignmentSeqLib::GetCigarData()
{
    SeqLib::Cigar a = alignment.GetCigar();
    cigar = CigarSeqLib(a);
    return cigar;
}

bool BamAlignmentSeqLib::HasSoftClips()
{
    return alignment.NumSoftClip()>0;
}

void BamAlignmentSeqLib::SetIsPaired(bool value)
{
    alignment.SetIsPaired(value);
}

void BamAlignmentSeqLib::SetIsProperPair(bool value)
{
    alignment.SetIsProperPair(value);
}

void BamAlignmentSeqLib::SetIsFirstMate(bool value)
{
    alignment.SetIsFirstMate(value);
}

void BamAlignmentSeqLib::SetIsMateMapped(bool value)
{
    alignment.SetIsMateMapped(value);
}

void BamAlignmentSeqLib::SetIsReverseStrand(bool value)
{
    alignment.SetIsReverseStrand(value);
}

void BamAlignmentSeqLib::SetIsMateReverseStrand(bool value)
{
    alignment.SetMateReverseFlag(value);
}

void BamAlignmentSeqLib::SetIsSecondMate(bool value)
{
    alignment.SetIsSecondMate(value);
}

void BamAlignmentSeqLib::SetQueryBases(std::string value)
{
    alignment.SetSequence(value);
    queryBases=value;
}

void BamAlignmentSeqLib::SetQualities(std::string value)
{
    alignment.SetQualities(value,0);
    qualities = value;
}

void BamAlignmentSeqLib::SetLength(int value)
{
    alignment.SetLength(value);
}

void BamAlignmentSeqLib::SetMateRefID(int value)
{
    alignment.SetMateChrID(value);
}

void BamAlignmentSeqLib::SetMatePosition(int value)
{
    alignment.SetMatePosition(value);
}

void BamAlignmentSeqLib::SetInsertSize(int value)
{
    alignment.SetInsertSize(value);
}

bool BamAlignmentSeqLib::HasTag(std::string & val)
{
    std::string s;
    return alignment.GetTag(val,s);
}

bool BamAlignmentSeqLib::AddZTag(std::string & tag, std::string & value)
{
    alignment.AddZTag(tag,value);
    return true;
}

bool BamAlignmentSeqLib::AddFloatTag(std::string & tag, float value)
{
    alignment.AddFloatTag(tag, value);
    return true;
}

bool BamAlignmentSeqLib::EditZTag(std::string & tag, std::string & value)
{
    alignment.RemoveTag(tag.c_str());
    alignment.AddZTag(tag, value);
    return true;
}

bool BamAlignmentSeqLib::EditFloatTag(std::string & tag, float value)
{
    alignment.RemoveTag(tag.c_str());
    alignment.AddFloatTag(tag, value);
    return true;
}

bool BamAlignmentSeqLib::GetZTag(std::string & tag, std::string & value)
{
    return alignment.GetZTag(tag,value);
}

bool BamAlignmentSeqLib::GetFloatTag(std::string & tag, float & value)
{
    return alignment.GetFloatTag(tag,value);
}

SeqLib::BamRecord &BamAlignmentSeqLib::GetSeqLibRecord()
{
    return alignment;
}
