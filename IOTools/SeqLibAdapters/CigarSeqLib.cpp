#include "CigarSeqLib.h"


CigarSeqLib::CigarSeqLib(){}

CigarSeqLib::CigarSeqLib(SeqLib::Cigar & cigar)
{
    this->cigar = cigar;
}

CigarIter& CigarSeqLib::createIterator()
{
    SeqLib::Cigar::iterator it = cigar.begin();
    SeqLib::Cigar::iterator end = cigar.end();
    iterator = CigarIterSeqLib(it,end);
    return iterator;
}

void CigarSeqLib::push_back(CigarOp & field)
{
    CigarOpSeqLib child = dynamic_cast<CigarOpSeqLib&>(field);
    cigar.add(child.field);
}

void CigarSeqLib::clear()
{
    cigar.Clear();
}



CigarIterSeqLib::CigarIterSeqLib(){}

CigarIterSeqLib::CigarIterSeqLib(SeqLib::Cigar::iterator & iterator,SeqLib::Cigar::iterator & end)
{
    this->iterator = iterator;
    this->end = end;
}

CigarOp& CigarIterSeqLib::Next()
{
    SeqLib::CigarField field = *iterator;
    currentItem = CigarOpSeqLib(field);
    iterator++;
    return currentItem;
}

bool CigarIterSeqLib::HasNext()
{
    return iterator!=end;
}


CigarOpSeqLib::CigarOpSeqLib(){}

CigarOpSeqLib::CigarOpSeqLib(SeqLib::CigarField & field)
{
    this->field = field;
}

CigarOpSeqLib::CigarOpSeqLib(char type, uint32_t length)
{
    field = SeqLib::CigarField(type,length);
}

char CigarOpSeqLib::GetType()
{
    return field.Type();
}

int CigarOpSeqLib::GetLength()
{
    return (int) field.Length();
}

void CigarOpSeqLib::setType(char newType)
{
    field = SeqLib::CigarField(newType, field.Length());
}
