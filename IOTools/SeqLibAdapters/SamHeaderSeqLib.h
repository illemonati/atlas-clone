#ifndef SAMHEADERSEQLIB_H
#define SAMHEADERSEQLIB_H

#include <vector>

#include "SeqLib/BamHeader.h"

#include "../IOAbstractClasses/SamHeader.h"


//Sequence
class SamSequenceSeqLib : public SamSequence{
private:
    std::string name;
    int length;
public:
    SamSequenceSeqLib();
    SamSequenceSeqLib(std::string&, int&);

    std::string GetName();
    int GetLength();

    ~SamSequenceSeqLib(){}
};

class SamSequenceIterSeqLib : public SamSequenceIter{
private:
    std::vector<SamSequenceSeqLib>::iterator iterator;
    std::vector<SamSequenceSeqLib>::iterator end;
    SamSequenceSeqLib currentItem;
public:
    SamSequenceIterSeqLib();
    SamSequenceIterSeqLib(std::vector<SamSequenceSeqLib>::iterator,std::vector<SamSequenceSeqLib>::iterator);

    SamSequence& Next();
    bool HasNext();

    ~SamSequenceIterSeqLib(){}
};

class SamSequenceDictionarySeqLib : public SamSequenceDictionary{
private:
    std::vector<SamSequenceSeqLib> sequences;
    SamSequenceIterSeqLib iterator;
public:
    SamSequenceDictionarySeqLib();
    SamSequenceDictionarySeqLib(std::vector<SamSequenceSeqLib>&);

    bool Contains(std::string&);
    int Size();
    SamSequenceIter* CreateIterator();
    SamSequenceIter *CreateIteratorEnd(int);

    ~SamSequenceDictionarySeqLib(){}
};

//ReadGroup
class SamReadGroupSeqLib : public SamReadGroup{
private:
    std::string library;
    std::string platformUnit;
    std::string predictedInsertSize;
    std::string productionDate;
    std::string program;
    std::string sample;
    std::string sequencingCenter;
    std::string sequencingTechnology;
    std::string id;
public:
    SamReadGroupSeqLib();
    SamReadGroupSeqLib(std::string,std::string="",std::string="",std::string="",std::string="",std::string="",std::string="",std::string="",std::string="");

    //Getter
    std::string GetLibrary();
    std::string GetPlatformUnit();
    std::string GetPredictedInsertSize();
    std::string GetProductionDate();
    std::string GetProgram();
    std::string GetSample();
    std::string GetSequencingCenter();
    std::string GetSequencingTechnology();
    std::string GetID();

    //Setter
    void SetLibrary(std::string);
    void SetPlatformUnit(std::string);
    void SetPredictedInsertSize(std::string);
    void SetProductionDate(std::string);
    void SetProgram(std::string);
    void SetSample(std::string);
    void SetSequencingCenter(std::string);
    void SetSequencingTechnology(std::string);

    ~SamReadGroupSeqLib(){}
};

class SamReadGroupIterSeqLib : public SamReadGroupIter{
private:
    std::vector<SamReadGroupSeqLib>::iterator iterator;
    std::vector<SamReadGroupSeqLib>::iterator end;
    SamReadGroupSeqLib currentItem;
public:
    SamReadGroupIterSeqLib();
    SamReadGroupIterSeqLib(std::vector<SamReadGroupSeqLib>::iterator, std::vector<SamReadGroupSeqLib>::iterator);

    SamReadGroup &Next();
    SamReadGroup& Get(int&);
    bool HasNext();

    ~SamReadGroupIterSeqLib(){}
};

class SamReadGroupDictionarySeqLib : public SamReadGroupDictionary{
private:
    std::vector<SamReadGroupSeqLib> readGroups;
    SamReadGroupIterSeqLib iterator;
public:
    SamReadGroupDictionarySeqLib();
    SamReadGroupDictionarySeqLib(std::vector<SamReadGroupSeqLib>);

    bool Contains(std::string&);
    int Size();
    void Add(std::string&);
    void Clear();
    SamReadGroupIter& CreateIterator();

    ~SamReadGroupDictionarySeqLib(){}
};

//Header
class SamHeaderSeqLib : public SamHeader
{
private:
    SeqLib::BamHeader samHeader;
    SamSequenceDictionarySeqLib sequences;
    SamReadGroupDictionarySeqLib readGroups;
    std::string version;
    std::string groupOrder;
    std::string sortOrder;

     std::vector<std::string> explode(const std::string&, const char&); //Code from https://stackoverflow.com/questions/890164/how-can-i-split-a-string-by-a-delimiter-into-an-array
    void parseHeader(std::string);
    void fillHeader(std::vector<std::string>&);
    SamSequenceSeqLib fillSequence(std::vector<std::string>&);
    SamReadGroupSeqLib fillReadGroup(std::vector<std::string>&);

public:
    SamHeaderSeqLib();
    SamHeaderSeqLib(SamHeaderSeqLib& samHeader);
    SamHeaderSeqLib(SeqLib::BamHeader);

    //Getter
    std::string GetVersion();
    std::string GetGroupOrder();
    std::string GetSortOrder();

    SamSequenceDictionary& GetSequences();
    SamReadGroupDictionary& GetReadGroups();

    ~SamHeaderSeqLib(){}

    void SetSamHeader(SeqLib::BamHeader);
};

#endif // SAMHEADERSEQLIB_H
