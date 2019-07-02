#ifndef IOABAMREADER_H
#define IOABAMREADER_H

#include <string>
#include "SamHeader.h"
#include "BamAlignment.h"
#include "RefVector.h"

/** Stream in reads from BAM/SAM/(CRAM)  */
class BamReader
{
public:
    /**
     * Open a BAM/SAM/(CRAM) file for streaming
     * @param filepath to a BAM/SAM/(CRAM) file
     * @return true if open was successful
     */
    virtual bool Open(std::string& filepath)=0;
    /**
     * Close the streaming
     * @return true if close was successful
     */
    virtual bool Close()=0;
    /**
     * Looks in BAM file's directory for a matching index file
     * @return true if there is a matching index file
     */
    virtual bool LocateIndex()=0;
    /**
     * Get a copy of the header of the file
     * @return a pointer to a SamHeader
     */
    virtual SamHeader* GetHeader()=0;
    /**
     * Jump to a specific region
     * @param the sequence ID and the position in the sequence
     * @return true if the jump have been done successfully
     */
    virtual bool Jump(int refid, int pos)=0;
    /**
     * Retrieve the next read from the available input streams.
     * @param Read to fill with data
     * @return true if a next read is available
     */
    virtual bool GetNextAlignment(BamAlignment*& alignment)=0;
    /**
     * Returns the number of reference sequences
     * @return the number of reference sequences
     */
    virtual int GetReferenceCount()=0;
    /**
     * Returns all reference sequence entries
     * @return a RefVector
     */
    virtual RefVector& GetReferenceData()=0;
    /**
     * Creates an index file for current BAM file
     * @return true if the index has been successfully created
     */
    virtual bool CreateIndex()=0;
    /**
     * Rewind the read file
     */
    virtual void Rewind()=0;
    /**
     * Tell the position in the file
     * @return the position in the file
     */
    virtual int tell()=0;

    virtual ~BamReader(){}
};

#endif // IOABAMREADER_H
