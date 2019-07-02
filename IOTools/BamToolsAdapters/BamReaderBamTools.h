#ifndef BAMREADERBAMTOOLS_H
#define BAMREADERBAMTOOLS_H

#include <string>

#include "IOAbstractClasses/SamHeader.h"
#include "SamHeaderBamTools.h"
#include "IOAbstractClasses/BamAlignment.h"
#include "BamAlignmentBamTools.h"
#include "IOAbstractClasses/RefVector.h"
#include "RefVectorBamTools.h"
#include "IOAbstractClasses/BamReader.h"

#include "../bamtools/api/BamReader.h"

class BamReaderBamTools : public BamReader
{
private:
    BamTools::BamReader bamReader;
    RefVectorBamTools refVectors;
    BamAlignmentBamTools nextAlignment;
    SamHeaderBamTools samHeader;

public:
    BamReaderBamTools();

    bool Open(std::string& filepath);
    bool Close();
    bool LocateIndex();
    SamHeader* GetHeader();
    bool Jump(int refid, int pos);
    bool GetNextAlignment(BamAlignment*& alignment);
    int GetReferenceCount();
    RefVector& GetReferenceData();
    bool CreateIndex();
    void Rewind();
    int tell();

    ~BamReaderBamTools(){}
};

#endif // BAMREADERBAMTOOLS_H
