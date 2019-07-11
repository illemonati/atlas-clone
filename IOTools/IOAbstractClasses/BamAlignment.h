#ifndef IOABAMALIGNMENT_H
#define IOABAMALIGNMENT_H

#include <string>
#include <vector>
#include "Cigar.h"

/**
 * Class to store and interact with a SAM alignment record
 */

class BamAlignment
{
public:

    //Getter
    /** Get the chromosome ID of the read */
    virtual int GetRefID()=0;
    /** Get the mapping quality */
    virtual int GetMapQuality()=0;
    /** Get the number of query bases of this read */
    virtual int GetLength()=0;
    /** Get the alignment position */
    virtual int GetPosition()=0;
    /** Get the chromosome ID of the mate */
    virtual int GetMateRefID()=0;
    /** Get the alignment position of mate */
    virtual int GetMatePosition()=0;
    /** Get the insert size for this read */
    virtual int GetInsertSize()=0;
    /** 'aligned' sequence (QueryBases plus deletion, padding, clipping chars)*/
    virtual std::string GetAlignedBases()=0;
    /** 'original' sequence (contained in BAM file) */
    virtual std::string GetQueryBases()=0;
    /** Get the quality scores of this read as a string */
    virtual std::string GetQualities()=0;
    /** Get readable chromosome name */
    virtual std::string GetName()=0;
    /** Returns true if alignment mapped to reverse strand */
    virtual bool IsReverseStrand()=0;
    /** Returns true if alignment part of paired-end read */
    virtual bool IsPaired()=0;
    /** Returns true if alignment is second mate on read */
    virtual bool IsSecondMate()=0;
    /** Returns true if alignment is part of read that satisfied paired-end resolution */
    virtual bool IsProperPair()=0;
    /** Returns true if alignment is mapped */
    virtual bool IsMapped()=0;
    /** Returns true if this read failed quality control */
    virtual bool IsFailedQC()=0;
    /** Returns true if reported position is primary alignment */
    virtual bool IsPrimaryAlignment()=0;
    /** Returns true if alignment is supplementary */
    virtual bool IsSupplementary()=0;
    /** Returns true if this read is a PCR duplicate */
    virtual bool IsDuplicate()=0;
    /** Get CIGAR for this alignment */
    virtual Cigar& GetCigarData()=0;
    /** Returns true if the alignment contains soft clips */
    virtual bool HasSoftClips()=0;

    //Setter
    /** Set value of "alignment part of paired-end read" flag */
    virtual void SetIsPaired(bool)=0;
    /** Set value of "alignment is part of read that satisfied paired-end resolution" flag */
    virtual void SetIsProperPair(bool)=0;
    /** Set value of "alignment is first mate" flag */
    virtual void SetIsFirstMate(bool)=0;
    /** Set value of "alignment's mate is mapped" flag */
    virtual void SetIsMateMapped(bool)=0;
    /** Set value of "alignment mapped to reverse strand" flag */
    virtual void SetIsReverseStrand(bool)=0;
    /** Set value of "alignment's mate mapped to reverse strand" flag */
    virtual void SetIsMateReverseStrand(bool)=0;
    /** Set value of "alignment is second mate on read" flag */
    virtual void SetIsSecondMate(bool)=0;
    /** Set value of the query base */
    virtual void SetQueryBases(std::string)=0;
    /** Set value of the qualities */
    virtual void SetQualities(std::string)=0;
    /** Set value of the length of alignment */
    virtual void SetLength(int)=0;
    /** Set value of the mate ID */
    virtual void SetMateRefID(int)=0;
    /** Set value of the mate position */
    virtual void SetMatePosition(int)=0;
    /** Set value of the insert size */
    virtual void SetInsertSize(int)=0;

    //Tag
    /**
     * Check if the alignment contains a specfic tag
     * @param the name of the tag
     * @return true if the alignment contains the tag
     */
    virtual bool HasTag(std::string&)=0;
    /**
     * Add a string Tag (Z) to the alignment
     * @param the name of the tag and its value
     * @return true if the tag has been successfully added
     */
    virtual bool AddZTag(std::string&,std::string&)=0;
    /**
     * Add a float Tag (f) to the alignment
     * @param the name of the tag and its value
     * @return true if the tag has been successfully added
     */
    virtual bool AddFloatTag(std::string&, float)=0;
    /**
     * Change the value of a string Tag (Z) to the alignment
     * @param the name of the tag to edit and its new value
     * @return true if the tag has been successfully edited
     */
    virtual bool EditZTag(std::string&,std::string&)=0;
    /**
     * Change the value of a floag Tag (f) to the alignment
     * @param the name of the tag to edit and its new value
     * @return true if the tag has been successfully edited
     */
    virtual bool EditFloatTag(std::string&, float)=0;
    /**
     * Get the value of a string Tag (Z) of the alignment
     * @param the name of the tag and the reference value
     * @return true if the tag has been found in the alignment
     */
    virtual bool GetZTag(std::string&, std::string&)=0;
    /**
     * Get the value of a float Tag (f) of the alignment
     * @param the name of the tag and the reference value
     * @return true if the tag has been found in the alignment
     */
    virtual bool GetFloatTag(std::string&, float&)=0;

    virtual ~BamAlignment(){}
};

#endif // IOABAMALIGNMENT_H
