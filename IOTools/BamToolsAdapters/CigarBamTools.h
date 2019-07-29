#ifndef CIGARBAMTOOLS_H
#define CIGARBAMTOOLS_H

#include <vector>

#include "../IOAbstractClasses/Cigar.h"

#include "bamtools/api/BamAux.h"

class CigarOpBamTools : public CigarOp
{
private:
    BamTools::CigarOp field;
public:
    CigarOpBamTools();
    CigarOpBamTools(BamTools::CigarOp&);
    CigarOpBamTools(char,uint32_t);

    char GetType();
    int GetLength();

    void setType(char newType);

    ~CigarOpBamTools(){}
};

class CigarIterBamTools : public CigarIter
{
private:
    std::vector<BamTools::CigarOp>::iterator iterator;
    std::vector<BamTools::CigarOp>::iterator end;
    CigarOpBamTools currentItem;
public:
    CigarIterBamTools();
    CigarIterBamTools(std::vector<BamTools::CigarOp>::iterator& iterator, std::vector<BamTools::CigarOp>::iterator& end);

    CigarOp& Next();
    bool HasNext();

    ~CigarIterBamTools(){}
};

class CigarBamTools : public Cigar
{
private:
    std::vector<BamTools::CigarOp> cigar;
    CigarIterBamTools iterator;
public:
    CigarBamTools();
    CigarBamTools(std::vector<BamTools::CigarOp>&);

    CigarIter& createIterator();

    void push_back(CigarOp&);
    void clear();

    ~CigarBamTools(){}
};

#endif // CIGARBAMTOOLS_H
