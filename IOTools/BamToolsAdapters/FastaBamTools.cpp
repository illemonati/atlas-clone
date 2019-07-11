#include "FastaBamTools.h"



bool FastaBamTools::GetSequence(const int &refid, const int &start, const int &stop, std::string &ref)
{
    return reference.GetSequence(refid,start,stop,ref);
}

bool FastaBamTools::Open(std::string &filepath, std::string &index)
{
    return reference.Open(filepath,index);;
}

bool FastaBamTools::GetBase(const int &refid, const int &position, char &base)
{
    return reference.GetBase(refid,position,base);
}
