#include "CigarBamTools.h"

CigarOpBamTools::CigarOpBamTools(){}

CigarOpBamTools::CigarOpBamTools(BamTools::CigarOp & field)
{
    this->field = field;
}

CigarOpBamTools::CigarOpBamTools(char type, uint32_t length)
{
    field = BamTools::CigarOp(type,length);
}

char CigarOpBamTools::GetType()
{
    return field.Type;
}

int CigarOpBamTools::GetLength()
{
    return field.Length;
}



CigarIterBamTools::CigarIterBamTools(){}

CigarIterBamTools::CigarIterBamTools(std::vector<BamTools::CigarOp>::iterator &iterator,std::vector<BamTools::CigarOp>::iterator& end){
    this->iterator = iterator;
    this->end= end;
}


CigarOp& CigarIterBamTools::Next()
{
    BamTools::CigarOp op = *iterator;
    currentItem= CigarOpBamTools(op);
    iterator++;
    return currentItem;
}

bool CigarIterBamTools::HasNext()
{
    return iterator!=end;
}



CigarBamTools::CigarBamTools(){}

CigarBamTools::CigarBamTools(std::vector<BamTools::CigarOp> & cigar)
{
    this->cigar=cigar;
}

CigarIter &CigarBamTools::createIterator()
{
    std::vector<BamTools::CigarOp>::iterator it = cigar.begin();
    std::vector<BamTools::CigarOp>::iterator end = cigar.end();
    iterator = CigarIterBamTools(it,end);
    return iterator;
}

void CigarBamTools::push_back(CigarOp &)
{

}

void CigarBamTools::clear()
{

}
