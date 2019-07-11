#include "IOTool.h"

IOTool* IOTool::instance = nullptr;
std::string IOTool::tool = "";

IOTool* IOTool::getInstance(){
    if(!instance){
        if(tool=="seqlib"){
#ifndef SEQLIB
            std::cout << "Unable to SeqLib because the program has been compiled without it.\n";
            exit(0);
#else
            std::cout << "Using SeqLib as Writing/Reading Tool." << std::endl;
            instance=new SeqLibFactory();
#endif
        }else{
            std::cout << "Using BamTools as Writing/Reading Tool." << std::endl;
            instance=new BamToolsFactory();
        }
    }
    return instance;
}

void IOTool::setParameter(std::string param){
    tool = param;
    instance = nullptr;
}

//BamToolsFactory
BamAlignment* BamToolsFactory::createBamAlignment(){
    return new BamAlignmentBamTools;
}

BamReader* BamToolsFactory::createBamReader(){
    int i = 0;
    return new BamReaderBamTools;
}

BamWriter* BamToolsFactory::createBamWriter(){
    return new BamWriterBamTools;
}

CigarOp* BamToolsFactory::createCigarOp(char type,uint32_t length){
    return new CigarOpBamTools(type, length);
}

Fasta* BamToolsFactory::createFasta(){
    return new FastaBamTools;
}

RefData* BamToolsFactory::createRefData(std::string& name, int length){
    return new RefDataBamTools(name, length);
}

RefVector* BamToolsFactory::createRefVector(){
    return new RefVectorBamTools;
}

SamHeader* BamToolsFactory::createSamHeader(SamHeader& samheader){
    SamHeaderBamTools child = dynamic_cast<SamHeaderBamTools&>(samheader);
    return new SamHeaderBamTools(child);
}

#ifdef SEQLIB
//SeqLibFactory
BamAlignment* SeqLibFactory::createBamAlignment(){
    return new BamAlignmentSeqLib;
}

BamReader* SeqLibFactory::createBamReader(){
    return new BamReaderSeqLib;
}

BamWriter* SeqLibFactory::createBamWriter(){
    return new BamWriterSeqLib;
}

CigarOp* SeqLibFactory::createCigarOp(char type,uint32_t length){
    return new CigarOpSeqLib(type, length);
}

Fasta* SeqLibFactory::createFasta(){
    return new FastaSeqLib;
}

RefData* SeqLibFactory::createRefData(std::string& name, int length){
    return new RefDataSeqLib(name, length);
}

RefVector* SeqLibFactory::createRefVector(){
    return new RefVectorSeqLib;
}

SamHeader* SeqLibFactory::createSamHeader(SamHeader& samheader){
    SamHeaderSeqLib child = dynamic_cast<SamHeaderSeqLib&>(samheader);
    return new SamHeaderSeqLib(child);
}
#endif
