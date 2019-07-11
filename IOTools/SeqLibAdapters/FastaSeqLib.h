#ifndef FASTASEQLIB_H
#define FASTASEQLIB_H


#include <string>
#include <vector>

#include "../IOAbstractClasses/Fasta.h"

#include "SeqLib/FastqReader.h"
#include "SeqLib/UnalignedSequence.h"

class FastaSeqLib : public Fasta
{
private:
    SeqLib::FastqReader reference;
    std::vector<std::string> buffer;
    int last_buffer_seq=-1;
public:
    bool GetSequence(const int& refid, const int& start, const int& stop, std::string& ref);
    bool Open(std::string& filepath, std::string& index);
    bool GetBase(const int& refid, const int& position, char& base);

    ~FastaSeqLib(){}
};

#endif // FASTASEQLIB_H
