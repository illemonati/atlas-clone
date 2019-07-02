#ifndef FASTABAMTOOLS_H
#define FASTABAMTOOLS_H

#include <string>

#include "../IOAbstractClasses/Fasta.h"

#include "bamtools/utils/bamtools_fasta.h"

class FastaBamTools : public Fasta
{
private:
    BamTools::Fasta reference;
public:
    bool GetSequence(const int& refid, const int& start, const int& stop, std::string& ref);
    bool Open(std::string& filepath, std::string& index);
    bool GetBase(const int& refid, const int& position, char& base);

    ~FastaBamTools(){}
};

#endif // FASTABAMTOOLS_H
