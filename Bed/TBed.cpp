#include "TBed.h"

//TBedChromosome
TBedChromosome::TBedChromosome(std::string &Name){
    name = Name;
}

TBedChromosome::~TBedChromosome(){
    //delete all windows
    windows.clear();
}

void TBedChromosome::addWindow(long start, long end){
    if(end <= start) throw "Window  [" + toString(start) + ", " + toString(end) + "] is not valid!";
    //check if there is overlap with any other window
    for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt){
        if((start >= windowIt->first && start < windowIt->second) || (end > windowIt->first && end <= windowIt->second))
            throw "Error reading bed file: window [" + toString(start) + ", " + toString(end) + ") overlaps window [" + toString(windowIt->first) + ", " + toString(windowIt->second) + ")!";
    }
    windows.insert(std::pair<long,long>(start, end));
}

void TBedChromosome::print(){
    std::cout << "Chromosome '" << name << "':" << std::endl;
    for(windowIt=windows.begin(); windowIt!=windows.end(); ++windowIt){
        std::cout << " - [" << windowIt->first << ", " << windowIt->second << "]" << std::endl;
    }
}

long TBedChromosome::size(){
    return windows.size();
}

bool TBedChromosome::begin(){
    windowIt = windows.begin();
    if(windowIt == windows.end()) return false;
    return true;
}

bool TBedChromosome::next(){
    ++windowIt;
    if(windowIt == windows.end()) return false;
    return true;
}

long TBedChromosome::curStart(){
    return windowIt->first;
}

long TBedChromosome::curEnd(){
    return windowIt->second;
}


//TBed
void TBed::readFile(){
    //open file
    std::istream* myStream = nullptr;
    if(filename.find(".gz")) myStream = new gz::igzstream(filename.c_str());
    else myStream = new std::ifstream(filename.c_str());
    if(!*myStream) throw "Failed to open BED file '" + filename + "'!";

    //tmp variables
    long lineNum = 0;
    std::vector<std::string> vec;
    curChr = "";

    //read file
    while(myStream->good() && !myStream->eof()){
        ++lineNum;
        std::string line;
        std::getline(*myStream, line);
        fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
        //skip empty lines
        if(vec.size() > 0){
            if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";

            //get chromosome
            if(vec[0] != curChr){
                chrIt = chromosomes.find(vec[0]);
                if(chrIt == chromosomes.end()){
                    chromosomes.insert(std::pair<std::string, TBedChromosome*>(vec[0], new TBedChromosome(vec[0])));
                    chrIt = chromosomes.find(vec[0]);
                }
                curChr = vec[0];
            }

            //add positions
            chrIt->second->addWindow(stringToLong(vec[1]), stringToLong(vec[2]));
        }
    }

    //close file
    delete myStream;
}

TBed::TBed(std::string Filename){
    filename = Filename;
    readFile();
    curChr = "";
}

TBed::~TBed(){
    //delete all chromosomes
    for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt){
        delete chrIt->second;
    }
    chromosomes.clear();
}

void TBed::setChr(const std::string &chr){
    curChr = chr;
    chrIt = chromosomes.find(curChr);
    if(chrIt!=chromosomes.end())
        chrIt->second->begin();
}

bool TBed::nextWindow(){
    if(chrIt==chromosomes.end()) return false;
    return chrIt->second->next();
}

long TBed::curWindowStart(){
    if(chrIt==chromosomes.end()) return -1;
    return chrIt->second->curStart();
}

long TBed::curWindowEnd(){
    if(chrIt==chromosomes.end()) return -1;
    return chrIt->second->curEnd();
}

void TBed::print(){
    std::cout << "Bed File '" << filename << "':" << std::endl;
    for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
        chrIt->second->print();
}

long TBed::size(){
    long s=0;
    for(chrIt=chromosomes.begin(); chrIt!=chromosomes.end(); ++chrIt)
        s += chrIt->second->size();
    return s;
}

int TBed::getNumChromosomes(){
    return chromosomes.size();
}

int TBed::getNumWindowsOnCurChr(){
    if(chrIt==chromosomes.end()) return 0;
    return chrIt->second->size();
}
