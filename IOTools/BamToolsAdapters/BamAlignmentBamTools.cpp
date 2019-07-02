#include "BamAlignmentBamTools.h"
#include <iostream>

BamAlignmentBamTools::BamAlignmentBamTools(){}

BamAlignmentBamTools::BamAlignmentBamTools(BamTools::BamAlignment alignment)
{
    this->alignment = alignment;
}

int BamAlignmentBamTools::GetRefID()
{
    return alignment.RefID;
}

int BamAlignmentBamTools::GetMapQuality()
{
    return alignment.MapQuality;
}

int BamAlignmentBamTools::GetLength()
{
    return alignment.Length;
}

int BamAlignmentBamTools::GetPosition()
{
    return alignment.Position;
}

int BamAlignmentBamTools::GetMateRefID()
{
    return alignment.MateRefID;
}

int BamAlignmentBamTools::GetMatePosition()
{
    return alignment.MatePosition;
}

int BamAlignmentBamTools::GetInsertSize()
{
    return alignment.InsertSize;
}

std::string BamAlignmentBamTools::GetAlignedBases()
{
    return alignment.AlignedBases;
}

std::string BamAlignmentBamTools::GetQueryBases()
{
    return alignment.QueryBases;
}

std::string BamAlignmentBamTools::GetQualities()
{
    return alignment.Qualities;
}

std::string BamAlignmentBamTools::GetName()
{
    return alignment.Name;
}

bool BamAlignmentBamTools::IsReverseStrand()
{
    return alignment.IsReverseStrand();
}

bool BamAlignmentBamTools::IsPaired()
{
    return alignment.IsPaired();
}

bool BamAlignmentBamTools::IsSecondMate()
{
    return alignment.IsSecondMate();
}

bool BamAlignmentBamTools::IsProperPair()
{
    return alignment.IsProperPair();
}

bool BamAlignmentBamTools::IsMapped()
{
    return alignment.IsMapped();
}

bool BamAlignmentBamTools::IsFailedQC()
{
    return alignment.IsFailedQC();
}

bool BamAlignmentBamTools::IsPrimaryAlignment()
{
    return alignment.IsPrimaryAlignment();
}

bool BamAlignmentBamTools::IsSupplementary()
{
    return alignment.IsSupplementary();
}

bool BamAlignmentBamTools::IsDuplicate()
{
    return alignment.IsDuplicate();
}

Cigar& BamAlignmentBamTools::GetCigarData()
{
    cigar = CigarBamTools(alignment.CigarData);
    return cigar;
}

bool BamAlignmentBamTools::HasSoftClips()
{
    std::vector<int> clipSizes;
    std::vector<int> readPositions;
    std::vector<int> genomePositions;
    return alignment.GetSoftClips(clipSizes,readPositions,genomePositions);
}

void BamAlignmentBamTools::SetIsPaired(bool value)
{
    alignment.SetIsPaired(value);
}

void BamAlignmentBamTools::SetIsProperPair(bool value)
{
    alignment.SetIsProperPair(value);
}

void BamAlignmentBamTools::SetIsFirstMate(bool value)
{
    alignment.SetIsFirstMate(value);
}

void BamAlignmentBamTools::SetIsReverseStrand(bool value)
{
    alignment.SetIsReverseStrand(value);
}

void BamAlignmentBamTools::SetIsMateMapped(bool value)
{
    alignment.SetIsMateMapped(value);
}

void BamAlignmentBamTools::SetIsMateReverseStrand(bool value)
{
    alignment.SetIsMateReverseStrand(value);
}

void BamAlignmentBamTools::SetIsSecondMate(bool value)
{
    alignment.SetIsSecondMate(value);
}

void BamAlignmentBamTools::SetQueryBases(std::string value)
{
    alignment.QueryBases = value;
}

void BamAlignmentBamTools::SetQualities(std::string value)
{
    alignment.Qualities = value;
}

void BamAlignmentBamTools::SetLength(int value)
{
    alignment.Length = value;
}

void BamAlignmentBamTools::SetMateRefID(int value)
{
    alignment.MateRefID = value;
}

void BamAlignmentBamTools::SetMatePosition(int value)
{
    alignment.MatePosition = value;
}

void BamAlignmentBamTools::SetInsertSize(int value)
{
    alignment.InsertSize = value;
}

bool BamAlignmentBamTools::HasTag(std::string & value)
{
    return alignment.HasTag(value);
}

bool BamAlignmentBamTools::AddZTag(std::string & tag, std::string & value)
{
    return alignment.AddTag(tag,"Z",value);
}

bool BamAlignmentBamTools::AddFloatTag(std::string & tag, float value)
{
    return alignment.AddTag(tag, "f", value);
}

bool BamAlignmentBamTools::EditZTag(std::string & tag, std::string & value)
{
    return alignment.EditTag(tag,"Z",value);
}

bool BamAlignmentBamTools::EditFloatTag(std::string & tag, float value)
{
    return alignment.EditTag(tag,"f",value);
}

bool BamAlignmentBamTools::GetZTag(std::string & tag, std::string & dest)
{
    return alignment.GetTag(tag, dest);
}

bool BamAlignmentBamTools::GetFloatTag(std::string & tag, float & dest)
{
    return alignment.GetTag(tag, dest);
}

BamTools::BamAlignment& BamAlignmentBamTools::GetBamToolsAlignment()
{
    return alignment;
}

