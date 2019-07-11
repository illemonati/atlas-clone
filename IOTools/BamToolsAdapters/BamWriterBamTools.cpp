#include "BamWriterBamTools.h"


bool BamWriterBamTools::Open(std::string &filepath, SamHeader &header, RefVector &references)
{
    RefVectorBamTools refVector = dynamic_cast<RefVectorBamTools&>(references);
    SamHeaderBamTools samHeader = dynamic_cast<SamHeaderBamTools&>(header);
    return bamWriter.Open(filepath,samHeader.GetSamHeader(),refVector.GetRefVector());
}

bool BamWriterBamTools::Close()
{
    bamWriter.Close();
    return true;
}

bool BamWriterBamTools::SaveAlignment(BamAlignment &alignment)
{
    auto & child = dynamic_cast<BamAlignmentBamTools&>(alignment);
    return bamWriter.SaveAlignment(child.GetBamToolsAlignment());
}
