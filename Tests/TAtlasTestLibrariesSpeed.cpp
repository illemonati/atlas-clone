#include "TAtlasTestLibrariesSpeed.h"

void TAtlasTestLibrarieSpeed::testSpeedBamReading()
{
    BamAlignment* alignment;

    for (int i=0; i < librariesName.size(); i++) {
        BamReader* reader = bamReaders.at(i);
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        reader->Open(bamFile);

        while(reader->GetNextAlignment(alignment)){
            alignment->GetQualities();
            alignment->GetQueryBases();
            alignment->GetAlignedBases();
            alignment->GetCigarData();
        }

        reader->Close();

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
        logfile->list("Tool["+librariesName.at(i)+"] time for reading all alignments : "+convertNanoSecondToReadableString(duration));
    }
}

void TAtlasTestLibrarieSpeed::testSpeedFastaReading()
{
    std::string indexes = fastaFile+".fai";

    //We first need to know the length of every sequences
    std::vector<long> length;
    BamReader* bamReader = bamReaders.at(0);
    bamReader->Open(bamFile);
    SamSequenceIter* iterator = bamReader->GetHeader()->GetSequences().CreateIterator();
    while(iterator->HasNext()){
        SamSequence& seq = iterator->Next();
        length.push_back(seq.GetLength());
    }

    //Then we run the speed test with a window of 1'000'000 bases
    logfile->list("Will request sequences for "+std::to_string(length.size())+" chromosome.");
    std::string sequence;
    char base;
    for (int i=0; i < librariesName.size(); i++) {
        Fasta* reader = fastaReaders.at(i);
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        reader->Open(fastaFile,indexes);


        for(unsigned i=0; i<length.size();i++){
            long nb_bases=length.at(i);
            if(nb_bases>0){
                long start=0;
                long end=1000000;
                while(end<nb_bases){
                    reader->GetSequence(i,start,end,sequence);
                    reader->GetBase(i,start,base);
                    start=end+1;
                    end+=1000000;
                }
                end=nb_bases;
                if(end>start){
                    reader->GetSequence(i,start,end,sequence);
                    reader->GetBase(i,start,base);
                }
            }
        }

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
        logfile->list("Tool["+librariesName.at(i)+"] time for reading all alignments : "+convertNanoSecondToReadableString(duration));
    }
}

std::string TAtlasTestLibrarieSpeed::convertNanoSecondToReadableString(long ms)
{
    int milliseconds = ms%1000;
    int seconds=(ms/1000)%60;
    int minutes=(ms/(1000*60))%60;
    std::string readableTime = std::to_string(minutes)+"m "+std::to_string(seconds)+"s "+std::to_string(milliseconds)+"ms";
    return readableTime;
}




TAtlasTestLibrarieSpeed::TAtlasTestLibrarieSpeed(TParameters &params, TLog *logfile):TAtlasTest(params, logfile){
    this->logfile=logfile;
    _name = "Librairie Compare Functions";
    bamFile=params.getParameterStringWithDefault("bam","");
    fastaFile=params.getParameterStringWithDefault("fasta","");
    librariesName = {"seqlib","bamtools"};
    for (std::vector<std::string>::iterator it = librariesName.begin(); it != librariesName.end(); ++it) {
        IOTool::setParameter(*it);
        bamReaders.push_back(IOTool::getInstance()->createBamReader());
        fastaReaders.push_back(IOTool::getInstance()->createFasta());
    }
}

bool TAtlasTestLibrarieSpeed::run()
{
    if(!bamFile.empty()){
        logfile->startIndent("Will run speed test for BAM Reading!");
        testSpeedBamReading();
        logfile->endIndent();
        if(!fastaFile.empty()){
            logfile->startIndent("Will run speed test for Fasta Reading!");
            testSpeedFastaReading();
            logfile->endIndent();
        }else{
            logfile->list("Couldn't run speed test for Fasta Reading because no file have been given in the parameters!");
        }
    }else{
        logfile->list("Couldn't run speed test for BAM Reading because no file have been given in the parameters! Won't run speed test for Fasta Reading because it needs SAM headers!");
    }

    for (auto p : bamReaders){
        delete p;
    }
    for (auto p : fastaReaders){
        delete p;
    }

    return true;
}

