#include "RefVectorBamTools.h"

RefVectorBamTools::RefVectorBamTools(BamTools::RefVector references)
{
    this->references=references;
}

RefVectorBamTools::RefVectorBamTools(){}

void RefVectorBamTools::push_back(RefData & refdata)
{
    RefDataBamTools child= dynamic_cast<RefDataBamTools&>(refdata);
    references.push_back(child.reference);
}

void RefVectorBamTools::clear()
{
    references.clear();
}

BamTools::RefVector RefVectorBamTools::GetRefVector()
{
    return references;
}



RefDataBamTools::RefDataBamTools(std::string& name, int length)
{
    reference = BamTools::RefData(name,length);
}
