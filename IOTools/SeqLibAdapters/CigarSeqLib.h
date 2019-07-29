#ifndef CIGARSEQLIB_H
#define CIGARSEQLIB_H

#include "../IOAbstractClasses/Cigar.h"

#include "SeqLib/BamRecord.h"

class CigarOpSeqLib : public CigarOp
{
private:
    SeqLib::CigarField field=SeqLib::CigarField('X',0);
public:
    CigarOpSeqLib();
    CigarOpSeqLib(SeqLib::CigarField&);
    CigarOpSeqLib(char, uint32_t);

    char GetType();
    int GetLength();
    void setType(char newType);

    ~CigarOpSeqLib(){}

    friend class CigarSeqLib;
};

class CigarIterSeqLib : public CigarIter
{
private:
    SeqLib::Cigar::iterator iterator;
    SeqLib::Cigar::iterator end;
    CigarOpSeqLib currentItem;
public:
    CigarIterSeqLib();
    CigarIterSeqLib(SeqLib::Cigar::iterator&,SeqLib::Cigar::iterator&);

    CigarOp& Next();
    bool HasNext();

    ~CigarIterSeqLib(){}
};

class CigarSeqLib : public Cigar
{
private:
    SeqLib::Cigar cigar;
    CigarIterSeqLib iterator;
public:
    CigarSeqLib();
    CigarSeqLib(SeqLib::Cigar&);

    CigarIter& createIterator();

    void push_back(CigarOp&);
    void clear();

    ~CigarSeqLib(){}
};

#endif // CIGARSEQLIB_H
