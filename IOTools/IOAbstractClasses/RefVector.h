#ifndef IOAREFVECTOR_H
#define IOAREFVECTOR_H
/**
 * Store a reference chromosome and its length
 */
class RefData
{
public:
    virtual ~RefData(){}
};

/**
 * List of reference chromosome
 */
class RefVector
{
public:
    /**
     * Add a RefData
     * @param the RefData to add
     */
    virtual void push_back(RefData& entry)=0;
    /**
     * Clear the RefVector
     */
    virtual void clear()=0;

    virtual ~RefVector(){}
};

#endif // IOAREFVECTOR_H
