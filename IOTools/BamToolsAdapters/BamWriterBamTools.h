#ifndef BAMWRITERBAMTOOLS_H
#define BAMWRITERBAMTOOLS_H

#include <string>
#include "../IOAbstractClasses/BamWriter.h"
#include "RefVectorBamTools.h"
#include "SamHeaderBamTools.h"
#include "BamAlignmentBamTools.h"


#include "../bamtools/api/BamWriter.h"


class BamWriterBamTools : public BamWriter
{
private:
    BamTools::BamWriter bamWriter;
public:

    bool Open(std::string& filepath, SamHeader& header, RefVector& references);
    bool Close();
    bool SaveAlignment(BamAlignment& alignment);

    ~BamWriterBamTools(){}
};

#endif // BAMWRITERBAMTOOLS_H
