#ifndef IOTOOL_H
#define IOTOOL_H

#include <string>
#include <iostream>
#include "IOAbstractClasses/BamAlignment.h"
#include "IOAbstractClasses/BamReader.h"
#include "IOAbstractClasses/BamWriter.h"
#include "IOAbstractClasses/Cigar.h"
#include "IOAbstractClasses/Fasta.h"
#include "IOAbstractClasses/RefVector.h"
#include "IOAbstractClasses/SamHeader.h"


#include "BamToolsAdapters/BamAlignmentBamTools.h"
#include "BamToolsAdapters/BamReaderBamTools.h"
#include "BamToolsAdapters/BamWriterBamTools.h"
#include "BamToolsAdapters/CigarBamTools.h"
#include "BamToolsAdapters/FastaBamTools.h"
#include "BamToolsAdapters/RefVectorBamTools.h"
#include "BamToolsAdapters/SamHeaderBamTools.h"


#include "SeqLibAdapters/BamAlignmentSeqLib.h"
#include "SeqLibAdapters/BamReaderSeqLib.h"
#include "SeqLibAdapters/BamWriterSeqLib.h"
#include "SeqLibAdapters/CigarSeqLib.h"
#include "SeqLibAdapters/FastaSeqLib.h"
#include "SeqLibAdapters/RefVectorSeqLib.h"
#include "SeqLibAdapters/SamHeaderSeqLib.h"

/**
 * Abstract Factory - Singelton
 * Purpose: Create input/output objects depending of the choosen tool
 */

class IOTool{
private:
    static IOTool* instance;
    static std::string tool;
public:
    /**
     * Singelton IOTool
     * @return an instance of IOTool
     */
    static IOTool* getInstance();

    /**
     * Set the tool to use
     * @param param The tool to use to create input/output objects
     */
    static void setParameter(std::string param);

    /**
     * Createa a BamAlignment
     * @return A pointer to a BamAlignment
     */
    virtual BamAlignment* createBamAlignment()=0;

    /**
     * Create a BamReader
     * @return A pointer to a BamReader
     */
    virtual BamReader* createBamReader()=0;
    /**
     * Create a BamWriter
     * @return A pointer to a BamWriter
     */
    virtual BamWriter* createBamWriter()=0;
    /**
     * Create a CigarOp
     * @param Type (MIDNSHP=X) and length
     * @return A pointer to a CigarOp
     */
    virtual CigarOp* createCigarOp(char,uint32_t)=0;
    /**
     * Create a Fasta reader
     * @return A pointer to a Fasta
     */
    virtual Fasta* createFasta()=0;
    /**
     * Create a RefData
     * @return A pointer to a RefData
     */
    virtual RefData* createRefData(std::string&, int)=0;
    /**
     * Create a RefVector
     * @return A pointer to a RefVector
     */
    virtual RefVector* createRefVector()=0;
    /**
     * Create a copy of SamHeader
     * @param samheader to copy
     * @return A pointer to a new SamHeader
     */
    virtual SamHeader* createSamHeader(SamHeader& samheader)=0;

    virtual ~IOTool(){}
};

class BamToolsFactory : public IOTool{
public:
    BamAlignment* createBamAlignment();
    BamReader* createBamReader();
    BamWriter*  createBamWriter();
    CigarOp* createCigarOp(char,uint32_t);
    Fasta* createFasta();
    RefData* createRefData(std::string&, int);
    RefVector* createRefVector();
    SamHeader* createSamHeader(SamHeader& samheader);

    ~BamToolsFactory(){}
};

class SeqLibFactory : public IOTool{
public:
    BamAlignment* createBamAlignment();
    BamReader* createBamReader();
    BamWriter*  createBamWriter();
    CigarOp* createCigarOp(char,uint32_t);
    Fasta* createFasta();
    RefData* createRefData(std::string&, int);
    RefVector* createRefVector();
    SamHeader* createSamHeader(SamHeader& samheader);

    ~SeqLibFactory(){}
};

#endif // IOTOOL_H
