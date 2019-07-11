#include "TAtlasTestLibraries.h"


TAtlasTestLibrarieCompareFunctions::TAtlasTestLibrarieCompareFunctions(TParameters &params, TLog *logfile):TAtlasTest(params, logfile){
    this->logfile=logfile;
    _name = "Librairie Compare Functions";
    lib1 = params.getParameterStringWithDefault("lib1","");
    lib2 = params.getParameterStringWithDefault("lib2","");
    std::unordered_set<std::string> libAvailables = {"seqlib","bamtools"}; //ADD HERE IF YOU WANT TO ALLOW OTHER LIBRAIRIES
    if(libAvailables.find(lib1) == libAvailables.end()){
        throw "Could not find lib1=["+lib1+"] in the available libraries";
    }
    if(libAvailables.find(lib2) == libAvailables.end()){
        throw "Could not find lib2=["+lib2+"] in the available libraries";
    }


    if(lib1.empty()||lib2.empty()){
        throw "Need two librairies to run librairiesCompareFunction. Must give the arguments lib1 and lib2.";
    }
    bamFile=params.getParameterStringWithDefault("bam","");
    fastaFile=params.getParameterStringWithDefault("fasta","");
}

bool TAtlasTestLibrarieCompareFunctions::run()
{
    if(fastaFile.empty() && bamFile.empty()){
        throw "Neither of bam or fasta file have been set in the parameters of the test. Must give the arguments bam or fasta to run at least one part of the test.";
    }

    bool bamReadingSuccessful=true;
    bool fastaReadingSuccessful=true;
    if(!bamFile.empty()){
        bamReadingSuccessful=compareBamReading();
    }else{
        logfile->list("Will run the test without bam reading comparison.");
    }


    if(!fastaFile.empty()){
        fastaReadingSuccessful=compareFastaReading();
    }else{
        logfile->list("Will run the test without fasta reading comparison.");
    }

    return fastaReadingSuccessful && bamReadingSuccessful;
}

bool TAtlasTestLibrarieCompareFunctions::compareBamReading()
{
    logfile->list("Create Fasta for lib1=["+lib1+"].");
    IOTool::setParameter(lib1);
    BamReader* lib1BamReader = IOTool::getInstance()->createBamReader();
    logfile->list("Create Fasta for lib2=["+lib2+"].");
    IOTool::setParameter(lib2);
    BamReader* lib2BamReader = IOTool::getInstance()->createBamReader();

    //Open Test
    bool lib1Open = lib1BamReader->Open(bamFile);
    bool lib2Open = lib2BamReader->Open(bamFile);

    if(lib1Open==lib2Open){
        logfile->list("Function [BamReader.Open] gave the same result!");
    }else{
        writeMismatch("BamReader.Open",std::to_string(lib1Open),std::to_string(lib2Open));
        logfile->list("Could not go on with the test because of open file failure!");
        return false;
    }
    if(!lib1Open){
        logfile->list("Could not go on with the test because both librairie have open file failure!");
        return true;
    }

    //Parsing Fasta
    BamAlignment* lib1Alignment;
    BamAlignment* lib2Alignment;
    lib1BamReader->GetNextAlignment(lib1Alignment);
    lib2BamReader->GetNextAlignment(lib2Alignment);

    checkResultMatchStr("BamAlignment.GetQueryBases",lib1Alignment->GetQueryBases(),lib2Alignment->GetQueryBases());
    checkResultMatchStr("BamAlignment.GetAlignedBases",lib1Alignment->GetAlignedBases(),lib2Alignment->GetAlignedBases());

    //Header
    SamHeader* lib1Header=lib1BamReader->GetHeader();
    SamHeader* lib2Header=lib2BamReader->GetHeader();
    checkResultMatchStr("SamHeader.GetVersion",lib1Header->GetVersion(),lib2Header->GetVersion());
    checkResultMatchStr("SamHeader.GetGroupOrder",lib1Header->GetGroupOrder(),lib2Header->GetGroupOrder());
    checkResultMatchStr("SamHeader.GetSortOrder",lib1Header->GetSortOrder(),lib2Header->GetSortOrder());

    //Sequences
    checkResultMatchInt("Sequences.Size",lib1Header->GetSequences().Size(),lib2Header->GetSequences().Size());
    checkResultMatchInt("Sequence.Length",lib1Header->GetSequences().CreateIterator()->Next().GetLength(),lib2Header->GetSequences().CreateIterator()->Next().GetLength());
    checkResultMatchStr("Sequence.Name",lib1Header->GetSequences().CreateIterator()->Next().GetName(),lib2Header->GetSequences().CreateIterator()->Next().GetName());

    //ReadGroup
    checkResultMatchInt("Sequences.Size",lib1Header->GetReadGroups().Size(),lib2Header->GetReadGroups().Size());

    //Others Test
    checkResultMatchInt("BamReader.GetReferenceCount",lib1BamReader->GetReferenceCount(),lib2BamReader->GetReferenceCount());

    //Close Test
    checkResultMatchInt("BamReader.Close",lib1BamReader->Close(),lib2BamReader->Close());

    return true;
}

bool TAtlasTestLibrarieCompareFunctions::compareFastaReading()
{
    logfile->list("Create Fasta reader for lib1=["+lib1+"].");
    IOTool::setParameter(lib1);
    Fasta* lib1Fasta = IOTool::getInstance()->createFasta();
    logfile->list("Create Fasta reader for lib2=["+lib2+"].");
    IOTool::setParameter(lib2);
    Fasta* lib2Fasta= IOTool::getInstance()->createFasta();

    //Open Test
    std::string indexes = fastaFile+".fai";
    bool lib1Open = lib1Fasta->Open(fastaFile,indexes);
    bool lib2Open = lib2Fasta->Open(fastaFile,indexes);

    if(lib1Open==lib2Open){
        logfile->list("Function [Fasta.Open] gave the same result!");
    }else{
        writeMismatch("Fasta.Open",std::to_string(lib1Open),std::to_string(lib2Open));
        logfile->list("Could not go on with the test because of open file failure!");
        return false;
    }
    if(!lib1Open){
        logfile->list("Could not go on with the test because both librairie have open file failure!");
        return true;
    }

    //GetSequence
    std::string lib1sequence;
    std::string lib2sequence;
    checkResultMatchInt("Fasta.GetSequence|HasSequence?", lib1Fasta->GetSequence(0,0,10,lib1sequence), lib2Fasta->GetSequence(0,0,10,lib2sequence));
    checkResultMatchStr("Fasta.GetSequence|Sequence",lib1sequence,lib2sequence);

    //GetBases
    char lib1base;
    char lib2base;
    checkResultMatchInt("Fasta.GetSequence|HasBase?", lib1Fasta->GetBase(0,3,lib1base), lib2Fasta->GetBase(0,3,lib2base));
    checkResultMatchStr("Fasta.GetSequence|Base",std::to_string(lib1base),std::to_string(lib2base));

    return true;
}

void TAtlasTestLibrarieCompareFunctions::checkResultMatchStr(std::string function, std::string resLib1, std::string resLib2)
{
    if(resLib1==resLib2){
        logfile->list("Function ["+function+"] gave the same result! Res="+resLib1);
    }else{
        writeMismatch(function,resLib1,resLib2);
    }
}

void TAtlasTestLibrarieCompareFunctions::checkResultMatchInt(std::string function, int resLib1, int resLib2)
{
    if(resLib1==resLib2){
        logfile->list("Function ["+function+"] gave the same result! Res="+std::to_string(resLib1));
    }else{
        writeMismatch(function,std::to_string(resLib1),std::to_string(resLib2));
    }
}

void TAtlasTestLibrarieCompareFunctions::writeMismatch(std::string function, std::string resLib1, std::string resLib2)
{
    logfile->list("Mismatch result for the function ["+function+"]. Lib1["+lib1+"] result is ["+resLib1+"] and Lib2["+lib2+"] result is ["+resLib2+"]");
}

