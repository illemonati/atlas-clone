#ifndef IOACIGAR_H
#define IOACIGAR_H

/**
 * Container for a single cigar operation
 */
class CigarOp
{
public:
    /**
     * Return the cigar op type as a char
     * @return a char (one of MIDNSHPX)
     */
    virtual char GetType()=0;
    virtual void setType(char newType)=0;
    /**
     * Return the length of the cigar op
     * @return a int representing the length
     */
    virtual int GetLength()=0;

    virtual ~CigarOp(){}
};

/**
 * Iterator to iterate through a Cigar
 */
class CigarIter
{
public:
    /**
     * Get the next cigar op
     * @return a CigarOp
     */
    virtual CigarOp& Next()=0;
    /**
     * Check if there is a next cigar op
     * @return true if there is a next cigar op
     */
    virtual bool HasNext()=0;

    virtual ~CigarIter(){}
};

/**
 * Cigar op vector
 */
class Cigar
{
public:
    /**
     * Create an iterator for the cigar
     * @return a CigarIter
     */
    virtual CigarIter& createIterator()=0;

    /**
     * Add a cigar op to the cigar
     * @param the cigar op to add
     */
    virtual void push_back(CigarOp&)=0;
    /**
     * Clear the whole cigar
     */
    virtual void clear()=0;

    virtual ~Cigar(){}
};



#endif // IOACIGAR_H
