#include "SamHeaderBamTools.h"


SamHeaderBamTools::SamHeaderBamTools(){}

SamHeaderBamTools::SamHeaderBamTools(SamHeaderBamTools &samHeader)
{
    this->samHeader = samHeader.samHeader;
}

SamHeaderBamTools::SamHeaderBamTools(BamTools::SamHeader samHeader)
{
    this->samHeader=samHeader;
}

std::string SamHeaderBamTools::GetVersion()
{
    return samHeader.Version;
}

std::string SamHeaderBamTools::GetGroupOrder()
{
    return samHeader.GroupOrder;
}

std::string SamHeaderBamTools::GetSortOrder()
{
    return samHeader.SortOrder;
}

SamSequenceDictionary& SamHeaderBamTools::GetSequences()
{
    sequences= SamSequenceDictionaryBamTools(samHeader.Sequences);
    return sequences;
}

SamReadGroupDictionary& SamHeaderBamTools::GetReadGroups()
{
    readGroups = SamReadGroupDictionaryBamTools(samHeader.ReadGroups);
    return readGroups;
}

BamTools::SamHeader SamHeaderBamTools::GetSamHeader()
{
    return samHeader;
}

/////////////////
//SamSequences
/////////////////
// Dictionary
SamSequenceDictionaryBamTools::SamSequenceDictionaryBamTools(){}

SamSequenceDictionaryBamTools::SamSequenceDictionaryBamTools(BamTools::SamSequenceDictionary& sequences){
    this->sequences = sequences;
}

bool SamSequenceDictionaryBamTools::Contains(std::string & value)
{
    return sequences.Contains(value);
}

int SamSequenceDictionaryBamTools::Size()
{
    return sequences.Size();
}

SamSequenceIter * SamSequenceDictionaryBamTools::CreateIterator()
{
    iterator = SamSequenceIterBamTools(sequences.Begin(), sequences.End());
    return &iterator;
}

SamSequenceIter * SamSequenceDictionaryBamTools::CreateIteratorEnd(int sub)
{
    iterator = SamSequenceIterBamTools(sequences.End()-sub, sequences.End());
    return &iterator;
}

// Iterator

SamSequenceIterBamTools::SamSequenceIterBamTools(){}

SamSequenceIterBamTools::SamSequenceIterBamTools(BamTools::SamSequenceIterator iterator, BamTools::SamSequenceIterator end)
{
    this->iterator = iterator;
    this->end = end;
}

SamSequence &SamSequenceIterBamTools::Next()
{
    BamTools::SamSequence sequence = *iterator;
    currentItem= SamSequenceBamTools(sequence);
    iterator++;
    return currentItem;
}

bool SamSequenceIterBamTools::HasNext()
{
    return iterator!=end;
}

// Sequence

SamSequenceBamTools::SamSequenceBamTools(){}

SamSequenceBamTools::SamSequenceBamTools(BamTools::SamSequence & sequence)
{
    this->sequence = sequence;
}

SamSequenceBamTools::SamSequenceBamTools(std::string & name, int & length)
{
    sequence = BamTools::SamSequence(name, length);
}

std::string SamSequenceBamTools::GetName()
{
    return sequence.Name;
}

int SamSequenceBamTools::GetLength()
{
    return std::stoi(sequence.Length);
}


//////////////////
//SamReadGroups
//////////////////
// Dictionary
SamReadGroupDictionaryBamTools::SamReadGroupDictionaryBamTools(){}

SamReadGroupDictionaryBamTools::SamReadGroupDictionaryBamTools(BamTools::SamReadGroupDictionary readGroups)
{
    this->readGroups = readGroups;
}


bool SamReadGroupDictionaryBamTools::Contains(std::string & value)
{
    return readGroups.Contains(value);
}

int SamReadGroupDictionaryBamTools::Size()
{
   return readGroups.Size();
}

void SamReadGroupDictionaryBamTools::Add(std::string & readgroup)
{
    readGroups.Add(readgroup);
}

void SamReadGroupDictionaryBamTools::Clear()
{
    return readGroups.Clear();
}

SamReadGroupIter &SamReadGroupDictionaryBamTools::CreateIterator()
{
    iterator = SamReadGroupIterBamTools(readGroups.Begin(), readGroups.End());
    return iterator;
}

//Iterator
SamReadGroupIterBamTools::SamReadGroupIterBamTools(){}

SamReadGroupIterBamTools::SamReadGroupIterBamTools(BamTools::SamReadGroupIterator iterator, BamTools::SamReadGroupIterator end)
{
    this->iterator = iterator;
    this->end = end;
}

SamReadGroup &SamReadGroupIterBamTools::Next()
{
    BamTools::SamReadGroup readgroup = *iterator;
    currentItem= SamReadGroupBamTools(readgroup);
    iterator++;
    return currentItem;
}

SamReadGroup &SamReadGroupIterBamTools::Get(int & id)
{
    BamTools::SamReadGroup readgroup = *(iterator+id);
    currentItem= SamReadGroupBamTools(readgroup);
    return currentItem;
}

bool SamReadGroupIterBamTools::HasNext()
{
    return iterator!=end;
}

// ReadGroup

SamReadGroupBamTools::SamReadGroupBamTools(){}

SamReadGroupBamTools::SamReadGroupBamTools(BamTools::SamReadGroup & readGroup)
{
    this->readGroup=readGroup;
}

std::string SamReadGroupBamTools::GetLibrary()
{
    return readGroup.Library;
}

std::string SamReadGroupBamTools::GetPlatformUnit()
{
    return readGroup.PlatformUnit;
}

std::string SamReadGroupBamTools::GetPredictedInsertSize()
{
    return readGroup.PredictedInsertSize;
}

std::string SamReadGroupBamTools::GetProductionDate()
{
    return readGroup.ProductionDate;
}

std::string SamReadGroupBamTools::GetProgram()
{
    return readGroup.Program;
}

std::string SamReadGroupBamTools::GetSample()
{
    return readGroup.Sample;
}

std::string SamReadGroupBamTools::GetSequencingCenter()
{
    return readGroup.SequencingCenter;
}

std::string SamReadGroupBamTools::GetSequencingTechnology()
{
    return readGroup.SequencingTechnology;
}

std::string SamReadGroupBamTools::GetID()
{
    return readGroup.ID;
}

void SamReadGroupBamTools::SetLibrary(std::string value)
{
    readGroup.Library = value;
}

void SamReadGroupBamTools::SetPlatformUnit(std::string value)
{
    readGroup.PlatformUnit = value;
}

void SamReadGroupBamTools::SetPredictedInsertSize(std::string value)
{
    readGroup.PredictedInsertSize = value;
}

void SamReadGroupBamTools::SetProductionDate(std::string value)
{
    readGroup.ProductionDate = value;
}

void SamReadGroupBamTools::SetProgram(std::string value)
{
    readGroup.Program = value;
}

void SamReadGroupBamTools::SetSample(std::string value)
{
    readGroup.Sample = value;
}

void SamReadGroupBamTools::SetSequencingCenter(std::string value)
{
    readGroup.SequencingCenter = value;
}

void SamReadGroupBamTools::SetSequencingTechnology(std::string value)
{
    readGroup.SequencingTechnology = value;
}

BamTools::SamReadGroup SamReadGroupBamTools::GetReadGroup()
{
    return readGroup;
}
