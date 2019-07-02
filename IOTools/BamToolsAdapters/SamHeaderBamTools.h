#ifndef SAMHEADERBAMTOOLS_H
#define SAMHEADERBAMTOOLS_H


#include "../IOAbstractClasses/SamHeader.h"

#include "../bamtools/api/SamHeader.h"
#include "../bamtools/api/SamSequence.h"
#include "../bamtools/api/SamSequenceDictionary.h"
#include "../bamtools/api/SamReadGroup.h"
#include "../bamtools/api/SamReadGroupDictionary.h"

class SamSequenceBamTools : public SamSequence{
private:
    BamTools::SamSequence sequence;
public:
    SamSequenceBamTools();
    SamSequenceBamTools(BamTools::SamSequence&);
    SamSequenceBamTools(std::string&, int&);

    std::string GetName();
    int GetLength();

    ~SamSequenceBamTools(){}
};

class SamSequenceIterBamTools : public SamSequenceIter{
private:
    BamTools::SamSequenceIterator iterator;
    BamTools::SamSequenceIterator end;
    SamSequenceBamTools currentItem;
public:
    SamSequenceIterBamTools();
    SamSequenceIterBamTools(BamTools::SamSequenceIterator,BamTools::SamSequenceIterator);

    SamSequence& Next();
    bool HasNext();

    ~SamSequenceIterBamTools(){}
};

class SamSequenceDictionaryBamTools : public SamSequenceDictionary{
private:
    BamTools::SamSequenceDictionary sequences;
    SamSequenceIterBamTools iterator;
public:
    SamSequenceDictionaryBamTools();
    SamSequenceDictionaryBamTools(BamTools::SamSequenceDictionary&);

    bool Contains(std::string&);
    int Size();
    SamSequenceIter* CreateIterator();
    SamSequenceIter* CreateIteratorEnd(int);

    ~SamSequenceDictionaryBamTools(){}
};



class SamReadGroupBamTools : public  SamReadGroup{
private:
    BamTools::SamReadGroup readGroup;
public:
    SamReadGroupBamTools();
    SamReadGroupBamTools(BamTools::SamReadGroup&);

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

    BamTools::SamReadGroup GetReadGroup();

    ~SamReadGroupBamTools(){}
};

class SamReadGroupIterBamTools : public SamReadGroupIter{
private:
    BamTools::SamReadGroupIterator iterator;
    BamTools::SamReadGroupIterator end;
    SamReadGroupBamTools currentItem;
public:
    SamReadGroupIterBamTools();
    SamReadGroupIterBamTools(BamTools::SamReadGroupIterator, BamTools::SamReadGroupIterator);

    SamReadGroup& Next();
    SamReadGroup& Get(int&);
    bool HasNext();

    ~SamReadGroupIterBamTools(){}
};

class SamReadGroupDictionaryBamTools : public SamReadGroupDictionary{
private:
    BamTools::SamReadGroupDictionary readGroups;
    SamReadGroupIterBamTools iterator;
public:
    SamReadGroupDictionaryBamTools();
    SamReadGroupDictionaryBamTools(BamTools::SamReadGroupDictionary);

    bool Contains(std::string&);
    int Size();
    void Add(std::string&);
    void Clear();
    SamReadGroupIter& CreateIterator();

    ~SamReadGroupDictionaryBamTools(){}
};


class SamHeaderBamTools : public SamHeader
{
private:
    BamTools::SamHeader samHeader;
    SamSequenceDictionaryBamTools sequences;
    SamReadGroupDictionaryBamTools readGroups;
public:

    SamHeaderBamTools();
    SamHeaderBamTools(SamHeaderBamTools& samHeader);
    SamHeaderBamTools(BamTools::SamHeader);

    //Getter
    std::string GetVersion();
    std::string GetGroupOrder();
    std::string GetSortOrder();
    SamSequenceDictionary& GetSequences();
    SamReadGroupDictionary& GetReadGroups();


    BamTools::SamHeader GetSamHeader();

    ~SamHeaderBamTools(){}
};

#endif // SAMHEADERBAMTOOLS_H
