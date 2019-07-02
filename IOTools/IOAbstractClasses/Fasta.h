#ifndef IOAFASTA_H
#define IOAFASTA_H

#include <string>

/**
 * Reader for fasta file
 */
class Fasta
{
public:
    /**
     * Get a subsequence of the fasta file
     * @param refid the reference ID
     * @param start the begin of the subsequence
     * @param stop the end of the subsequence
     * @param ref the string to fill with the subsequence
     * @return true if the sequence has been found
     */
    virtual bool GetSequence(const int& refid, const int& start, const int& stop, std::string& ref)=0;
    /**
     * Open a fasta file
     * @param filepath of the fasta and of the corresponding fasta index
     * @return true if the file has been successfully opened
     */
    virtual bool Open(std::string& filepath, std::string& index)=0;
    /**
     * Give the base of a specific position in the fasta file
     * @param reference ID and position
     * @param char to fill with the base
     * @return true if the base has been found
     */
    virtual bool GetBase(const int& refid, const int& position, char& base)=0;

    virtual ~Fasta(){}
};

#endif // IOAFASTA_H
