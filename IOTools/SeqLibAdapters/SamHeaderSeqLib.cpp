#include "SamHeaderSeqLib.h"

//Header

std::vector<std::string> SamHeaderSeqLib::explode(const std::string& str, const char& ch) {
    std::string next;
    std::vector<std::string> result;

    // For each character in the string
    for (std::string::const_iterator it = str.begin(); it != str.end(); it++) {
        // If we've hit the terminal character
        if (*it == ch) {
            // If we have some characters accumulated
            if (!next.empty()) {
                // Add them to the result vector
                result.push_back(next);
                next.clear();
            }
        } else {
            // Accumulate the next character into the sequence
            next += *it;
        }
    }
    if (!next.empty())
         result.push_back(next);
    return result;
}

void SamHeaderSeqLib::parseHeader(std::string header)
{
    std::vector<std::string> lines = explode(header, '\n');
    std::vector<std::string> fields;
    std::vector<SamSequenceSeqLib> seqs;
    std::vector<SamReadGroupSeqLib> rgs;
    for (size_t i = 0; i < lines.size(); i++) {
        fields = explode(lines[i],'\t');
        std::string& first_tag = fields.at(0);
        if(first_tag=="@HD"){
            fillHeader(fields);
        }else if(first_tag=="@SQ"){
            seqs.push_back(fillSequence(fields));
        }else if(first_tag=="@RG"){
            rgs.push_back(fillReadGroup(fields));
        }else{
            //else it's the PG tag and we're at the end of the BAM Header
            break;
        }
    }
    sequences = SamSequenceDictionarySeqLib(seqs);
    readGroups = SamReadGroupDictionarySeqLib(rgs);
}

inline void SamHeaderSeqLib::fillHeader(std::vector<std::string>& fields)
{
    for (size_t i = 0; i < fields.size(); i++) {
        std::string& field = fields.at(i);
        std::string tag = field.substr(0,2);
        if(tag=="VN"){
            version=field.substr(3);
        }else if(tag=="SO"){
            sortOrder=field.substr(3);
        }else if(tag=="GO"){
            groupOrder=field.substr(3);
        }
    }
}

inline SamSequenceSeqLib SamHeaderSeqLib::fillSequence(std::vector<std::string>& fields)
{
    std::string name;
    int length;
    for (size_t i = 0; i < fields.size(); i++) {
        std::string& field = fields.at(i);
        std::string tag = field.substr(0,2);
        if(tag=="SN"){
            name=field.substr(3);
        }else if(tag=="LN"){
            length=std::atoi(field.substr(3).c_str());
        }
    }
    return SamSequenceSeqLib(name,length);
}

inline SamReadGroupSeqLib SamHeaderSeqLib::fillReadGroup(std::vector<std::string>& fields)
{
    std::string library,platformUnit,predictedInsertSize,productionDate,program,sample,sequencingCenter,sequencingTechnology,id;
    for (size_t i = 0; i < fields.size(); i++) {
        std::string& field = fields.at(i);
        std::string tag = field.substr(0,2);
        if(tag=="ID"){
            id=field.substr(3);
        }else if(tag=="CN"){
            sequencingCenter=field.substr(3);
        }else if(tag=="LB"){
            library=field.substr(3);
        }else if(tag=="PG"){
            program=field.substr(3);
        }else if(tag=="PI"){
            predictedInsertSize=field.substr(3);
        }else if(tag=="DT"){
            productionDate=field.substr(3);
        }else if(tag=="PU"){
            platformUnit=field.substr(3);
        }else if(tag=="PL"){
            sequencingTechnology=field.substr(3);
        }else if(tag=="SM"){
            sample=field.substr(3);
        }
    }
    return SamReadGroupSeqLib(id,library,platformUnit,predictedInsertSize,productionDate,program,sample,sequencingCenter,sequencingTechnology);
}

SamHeaderSeqLib::SamHeaderSeqLib(){}

SamHeaderSeqLib::SamHeaderSeqLib(SamHeaderSeqLib &samHeader)
{
    this->samHeader=samHeader.samHeader;
    this->version=samHeader.version;
    this->sortOrder=samHeader.sortOrder;
    this->groupOrder=samHeader.groupOrder;
    this->sequences=samHeader.sequences;
    this->readGroups=samHeader.readGroups;
}

SamHeaderSeqLib::SamHeaderSeqLib(SeqLib::BamHeader samHeader)
{
    this->samHeader=samHeader;
    parseHeader(samHeader.AsString());
}

std::string SamHeaderSeqLib::GetVersion()
{
    return version;
}

std::string SamHeaderSeqLib::GetGroupOrder()
{
    return groupOrder;
}

std::string SamHeaderSeqLib::GetSortOrder()
{
    return sortOrder;
}

SamSequenceDictionary &SamHeaderSeqLib::GetSequences()
{
    return sequences;
}

SamReadGroupDictionary &SamHeaderSeqLib::GetReadGroups()
{
    return readGroups;
}

void SamHeaderSeqLib::SetSamHeader(SeqLib::BamHeader samHeader)
{
    this->samHeader=samHeader;
    parseHeader(samHeader.AsString());
}


//////////////
//SamSequence
/////////////
// Dictionary
SamSequenceDictionarySeqLib::SamSequenceDictionarySeqLib(){}

SamSequenceDictionarySeqLib::SamSequenceDictionarySeqLib(std::vector<SamSequenceSeqLib> & sequences)
{
    this->sequences = sequences;
}

bool SamSequenceDictionarySeqLib::Contains(std::string & value)
{
    for(std::vector<SamSequenceSeqLib>::iterator it = sequences.begin(); it != sequences.end(); ++it) {
        if(it->GetName()==value){
            return true;
        }
    }
    return false;
}

int SamSequenceDictionarySeqLib::Size()
{
    return sequences.size();
}

SamSequenceIter *SamSequenceDictionarySeqLib::CreateIterator()
{
    iterator = SamSequenceIterSeqLib(sequences.begin(), sequences.end());
    return &iterator;
}

SamSequenceIter *SamSequenceDictionarySeqLib::CreateIteratorEnd(int sub)
{
    iterator = SamSequenceIterSeqLib(sequences.end()-sub, sequences.end());
    return &iterator;
}

// Iterator
SamSequenceIterSeqLib::SamSequenceIterSeqLib(){}

SamSequenceIterSeqLib::SamSequenceIterSeqLib(std::vector<SamSequenceSeqLib>::iterator iterator, std::vector<SamSequenceSeqLib>::iterator end)
{
    this->iterator = iterator;
    this->end = end;
}

SamSequence &SamSequenceIterSeqLib::Next()
{
    currentItem = *iterator;
    iterator++;
    return currentItem;
}

bool SamSequenceIterSeqLib::HasNext()
{
    return iterator!=end;
}

// Sequence
SamSequenceSeqLib::SamSequenceSeqLib(){}

SamSequenceSeqLib::SamSequenceSeqLib(std::string & name, int & length)
{
    this->name = name;
    this->length = length;
}

std::string SamSequenceSeqLib::GetName()
{
    return name;
}

int SamSequenceSeqLib::GetLength()
{
    return length;
}


//////////////
//SamReadGroup
/////////////
// Dictionnary

SamReadGroupDictionarySeqLib::SamReadGroupDictionarySeqLib(){}

SamReadGroupDictionarySeqLib::SamReadGroupDictionarySeqLib(std::vector<SamReadGroupSeqLib> readGroups)
{
    this->readGroups=readGroups;
}

bool SamReadGroupDictionarySeqLib::Contains(std::string & value)
{
    for(std::vector<SamReadGroupSeqLib>::iterator it = readGroups.begin(); it != readGroups.end(); ++it) {
        if(it->GetID()==value){
            return true;
        }
    }
    return false;
}

int SamReadGroupDictionarySeqLib::Size()
{
    return readGroups.size();
}

void SamReadGroupDictionarySeqLib::Add(std::string & id)
{
    if(readGroups.empty() || !Contains(id)){
        readGroups.push_back(SamReadGroupSeqLib(id));
    }
}

void SamReadGroupDictionarySeqLib::Clear()
{
    readGroups.clear();
}

SamReadGroupIter &SamReadGroupDictionarySeqLib::CreateIterator()
{
    iterator = SamReadGroupIterSeqLib(readGroups.begin(), readGroups.end());
    return iterator;
}


// Iterator
SamReadGroupIterSeqLib::SamReadGroupIterSeqLib(){}

SamReadGroupIterSeqLib::SamReadGroupIterSeqLib(std::vector<SamReadGroupSeqLib>::iterator iterator, std::vector<SamReadGroupSeqLib>::iterator end)
{
    this->iterator=iterator;
    this->end=end;
}

SamReadGroup &SamReadGroupIterSeqLib::Next()
{
    currentItem = *iterator;
    iterator++;
    return currentItem;
}

SamReadGroup &SamReadGroupIterSeqLib::Get(int & id)
{
    currentItem= *(iterator+id);
    return currentItem;
}

bool SamReadGroupIterSeqLib::HasNext()
{
    return iterator!=end;
}

// ReadGroup

SamReadGroupSeqLib::SamReadGroupSeqLib(){}

SamReadGroupSeqLib::SamReadGroupSeqLib(std::string id, std::string library, std::string platformUnit, std::string predictedInsertSize, std::string productionDate,
                                       std::string program, std::string sample, std::string sequencingCenter, std::string sequencingTechnology)
{
    this->id=id;
    this->library=library;
    this->platformUnit=platformUnit;
    this->predictedInsertSize=predictedInsertSize;
    this->productionDate=productionDate;
    this->program=program;
    this->sample=sample;
    this->sequencingCenter=sequencingCenter;
    this->sequencingTechnology=sequencingTechnology;
}

std::string SamReadGroupSeqLib::GetLibrary()
{
    return library;
}

std::string SamReadGroupSeqLib::GetPlatformUnit()
{
    return platformUnit;
}

std::string SamReadGroupSeqLib::GetPredictedInsertSize()
{
    return predictedInsertSize;
}

std::string SamReadGroupSeqLib::GetProductionDate()
{
    return productionDate;
}

std::string SamReadGroupSeqLib::GetProgram()
{
    return program;
}

std::string SamReadGroupSeqLib::GetSample()
{
    return sample;
}

std::string SamReadGroupSeqLib::GetSequencingCenter()
{
    return sequencingCenter;
}

std::string SamReadGroupSeqLib::GetSequencingTechnology()
{
    return sequencingTechnology;
}

std::string SamReadGroupSeqLib::GetID()
{
    return id;
}

void SamReadGroupSeqLib::SetLibrary(std::string value)
{
    id=value;
}

void SamReadGroupSeqLib::SetPlatformUnit(std::string value)
{
    platformUnit=value;
}

void SamReadGroupSeqLib::SetPredictedInsertSize(std::string value)
{
    predictedInsertSize=value;
}

void SamReadGroupSeqLib::SetProductionDate(std::string value)
{
    productionDate=value;
}

void SamReadGroupSeqLib::SetProgram(std::string value)
{
    program=value;
}

void SamReadGroupSeqLib::SetSample(std::string value)
{
    sample=value;
}

void SamReadGroupSeqLib::SetSequencingCenter(std::string value)
{
    sequencingCenter=value;
}

void SamReadGroupSeqLib::SetSequencingTechnology(std::string value)
{
    sequencingTechnology=value;
}
