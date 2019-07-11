#ifndef IOASAMHEADER_H
#define IOASAMHEADER_H

#include <string>


/**
 * Container of a sequence
 */
class SamSequence{
public:
    /** Get the name of the sequence */
    virtual std::string GetName()=0;
    /** Get the length of the sequence */
    virtual int GetLength()=0;

    virtual ~SamSequence(){}
};

/**
 * Iterator of SamSequenceDictionary
 */
class SamSequenceIter{
public:
    /**
     * Get the next sequence
     * @return next SamSequence
     */
    virtual SamSequence& Next()=0;
    /**
     * @brief HasNext
     * @return true if there is a next sequence
     */
    virtual bool HasNext()=0;

    virtual ~SamSequenceIter(){}
};

/**
 * Vector of SamSequence
 */
class SamSequenceDictionary{
public:
    /**
     * Check if the vector contains a specific sequence name
     * @return true if find the sequence name
     */
    virtual bool Contains(std::string&)=0;
    /**
     * Get the size of the vector
     * @return
     */
    virtual int Size()=0;
    /**
     * Create an iterator from the begining
     * @return a SamSequenceIter
     */
    virtual SamSequenceIter* CreateIterator()=0;
    /**
     * Create an iterator from the end
     * @return a SamSequenceIter
     */
    virtual SamSequenceIter* CreateIteratorEnd(int=0)=0;

    virtual ~SamSequenceDictionary(){}
};

/**
 * Container for a  read group with setters and getters
 */
class SamReadGroup{
public:
    //Getter
    virtual std::string GetLibrary()=0;
    virtual std::string GetPlatformUnit()=0;
    virtual std::string GetPredictedInsertSize()=0;
    virtual std::string GetProductionDate()=0;
    virtual std::string GetProgram()=0;
    virtual std::string GetSample()=0;
    virtual std::string GetSequencingCenter()=0;
    virtual std::string GetSequencingTechnology()=0;
    virtual std::string GetID()=0;

    //Setter
    virtual void SetLibrary(std::string)=0;
    virtual void SetPlatformUnit(std::string)=0;
    virtual void SetPredictedInsertSize(std::string)=0;
    virtual void SetProductionDate(std::string)=0;
    virtual void SetProgram(std::string)=0;
    virtual void SetSample(std::string)=0;
    virtual void SetSequencingCenter(std::string)=0;
    virtual void SetSequencingTechnology(std::string)=0;

    virtual ~SamReadGroup(){}
};
/**
 * Iterator of SamReadGroupDictionary
 */

class SamReadGroupIter{
public:
    /**
     * Get the next read group
     * @return next SamReadGroup
     */
    virtual SamReadGroup& Next()=0;
    /**
     * Get a read group at a specific location
     * @return a SamReadGroup
     */
    virtual SamReadGroup& Get(int&)=0;
    /**
     * Check if there is a next read group
     * @return true if there is a next read group
     */
    virtual bool HasNext()=0;

    virtual ~SamReadGroupIter(){}
};

/**
 * Vector of Read Group
 */
class SamReadGroupDictionary{
public:
    /**
     * Check if the vector contains a specific read group id
     * @return true if the vector contains the read group id
     */
    virtual bool Contains(std::string&)=0;
    /**
     * Size of the read group vector
     * @return the number of read group
     */
    virtual int Size()=0;
    /**
     * Add a read group
     * @param the read group id to add
     */
    virtual void Add(std::string&)=0;
    /**
     * Clear the read group vector
     */
    virtual void Clear()=0;
    /**
     * Create an iterator to parse the read group vector
     * @return a SamReadGroupIter
     */
    virtual SamReadGroupIter& CreateIterator()=0;

    virtual ~SamReadGroupDictionary(){}
};

/**
 * Container of the header of a BAM/SAM file
 */
class SamHeader
{
public:
    //Getter
    /** Get the version format */
    virtual std::string GetVersion()=0;
    /** Get the group order */
    virtual std::string GetGroupOrder()=0;
    /** Get the sort order */
    virtual std::string GetSortOrder()=0;
    /** Get the list of sequence */
    virtual SamSequenceDictionary& GetSequences()=0;
    /** Get the list of read groups */
    virtual SamReadGroupDictionary& GetReadGroups()=0;

    virtual ~SamHeader(){}
};




#endif // IOASAMHEADER_H
