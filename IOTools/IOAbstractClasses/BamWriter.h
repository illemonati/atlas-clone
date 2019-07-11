#ifndef IOABAMWRITER_H
#define IOABAMWRITER_H

#include <string>
#include "SamHeader.h"
#include "RefVector.h"
#include "BamAlignment.h"

/**
 * Writer for BAM/SAM/(CRAM) file
 */
class BamWriter
{
public:
    /**
     * Open a BAM file for writing
     * @param filepath the name of the file to wrtie
     * @param header the header of the file to wrtie
     * @param references the references to write
     * @return true if the openning was successful
     */
    virtual bool Open(std::string& filepath, SamHeader& header, RefVector& references)=0;
    /**
     * Close the file
     * @return true if the file has been succefully closed
     */
    virtual bool Close()=0;
    /**
     * Write an alignment to the output file
     * @param alignment the alignment to write
     * @return true if the write was successful
     */
    virtual bool SaveAlignment(BamAlignment& alignment)=0;

    virtual ~BamWriter(){}
};

#endif // IOABAMWRITER_H
