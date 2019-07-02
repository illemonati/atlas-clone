#ifndef REFVECTORBAMTOOLS_H
#define REFVECTORBAMTOOLS_H

#include <string>

#include "../bamtools/api/BamAux.h"

#include "IOAbstractClasses/RefVector.h"

class RefVectorBamTools: public RefVector
{
private:
    BamTools::RefVector references;
public:
    RefVectorBamTools(BamTools::RefVector);
    RefVectorBamTools();

    void push_back(RefData&);
    void clear();

    BamTools::RefVector GetRefVector();

    ~RefVectorBamTools(){}
};

class RefDataBamTools : public RefData
{
private:
    BamTools::RefData reference;
public:
    RefDataBamTools(std::string&, int);

    ~RefDataBamTools(){}

    friend class RefVectorBamTools;
};

#endif // REFVECTORBAMTOOLS_H
