#ifndef BAMALIGNMENTBAMTOOLS_H
#define BAMALIGNMENTBAMTOOLS_H

#include <string>
#include <vector>
#include "../IOAbstractClasses/Cigar.h"
#include "CigarBamTools.h"
#include "../IOAbstractClasses/BamAlignment.h"

#include "bamtools/api/BamAlignment.h"

class BamWriterBamTools;

class BamAlignmentBamTools : public BamAlignment
{
private:
    BamTools::BamAlignment alignment;
    CigarBamTools cigar;
public:

    BamAlignmentBamTools();
    BamAlignmentBamTools(BamTools::BamAlignment);

    //Getter
    int GetRefID();
    int GetMapQuality();
    int GetLength();
    int GetPosition();
    int GetMateRefID();
    int GetMatePosition();
    int GetInsertSize();
    std::string GetAlignedBases();
    std::string GetQueryBases();
    std::string GetQualities();
    std::string GetName();
    bool IsReverseStrand();
    bool IsPaired();
    bool IsSecondMate();
    bool IsProperPair();
    bool IsMapped();
    bool IsFailedQC();
    bool IsPrimaryAlignment();
    bool IsSupplementary();
    bool IsDuplicate();
    Cigar& GetCigarData();
    bool HasSoftClips();

    //Setter
    void SetIsPaired(bool);
    void SetIsProperPair(bool);
    void SetIsFirstMate(bool);
    void SetIsMateMapped(bool);
    void SetIsReverseStrand(bool);
    void SetIsMateReverseStrand(bool);
    void SetIsSecondMate(bool);
    void SetQueryBases(std::string);
    void SetQualities(std::string);
    void SetLength(int);
    void SetMateRefID(int);
    void SetMatePosition(int);
    void SetInsertSize(int);


    //Tag
    bool HasTag(std::string&);
    bool AddZTag(std::string&,std::string&);
    bool AddFloatTag(std::string&, float);
    bool EditZTag(std::string&,std::string&);
    bool EditFloatTag(std::string&, float);
    bool GetZTag(std::string&, std::string&);
    bool GetFloatTag(std::string&, float&);

    BamTools::BamAlignment& GetBamToolsAlignment();

    ~BamAlignmentBamTools(){}
};

#endif // BAMALIGNMENTBAMTOOLS_H
