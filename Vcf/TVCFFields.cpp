#include "TVCFFields.h"

//TVCFField
TVCFField::TVCFField(const std::string &Tag, const std::string &Description){
    tag = Tag;
    description = Description;
    accepted = false;
    used = false;
}

void TVCFField::writeHeader(std::string type, gz::ogzstream &vcf){
    vcf << "##" << type << "=<ID=" << tag << ',' << description << ">\n";
}


//TVCFFieldVector
TVCFFieldVector::TVCFFieldVector(std::string Type){
    _type = Type;
}

std::string TVCFFieldVector::type(){ return _type; }

void TVCFFieldVector::add(std::string tag, std::string description){
    availableFields.emplace_back(tag, description);
}

TVCFField *TVCFFieldVector::getFieldPointer(const std::string &tag){
    for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it){
        if(tag == it->tag) return &(*it);
    }
    throw "VCF " + _type + " field '" + tag + "' is unknown!";
}

bool TVCFFieldVector::fieldExists(const std::string &tag){
    for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it){
        if(it->tag == tag)
            return true;
    }
    return false;
}

void TVCFFieldVector::acceptField(const std::string &tag){
    getFieldPointer(tag)->accepted = true;
}

bool TVCFFieldVector::useField(std::string tag){
    TVCFField* pt = getFieldPointer(tag);
    if(pt->accepted){
        if(!pt->used){
            usedFields.push_back(pt);
            pt->used = true;
        }
        return true;
    }
    return false;
}

void TVCFFieldVector::clearUsed(){
    usedFields.clear();
    for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it)
        it->used = false;
}

int TVCFFieldVector::numUsed(){
    return usedFields.size();
}

std::string TVCFFieldVector::getListOfUsedFields(std::string delim){
    std::string out;
    for(TVCFField* it : usedFields){
        if(out.length() > 0) out += delim;
        out += it->tag;
    }
    return out;
}

void TVCFFieldVector::fillVectorWithTagsOfUsedFields(std::vector<std::string> &vec){
    vec.clear();
    for(std::vector<TVCFField*>::iterator it = usedFields.begin(); it != usedFields.end(); ++it)
        vec.push_back((*it)->tag);
}

void TVCFFieldVector::writeVCFHeader(gz::ogzstream &vcf){
    //info fields
    for(std::vector<TVCFField*>::iterator it = usedFields.begin(); it != usedFields.end(); ++it)
        (*it)->writeHeader(_type, vcf);
}

void TVCFFieldVector::reportUsedFields(TLog *logfile){
    if(numUsed() > 0)
        logfile->list("Will print these VCF " + _type + " fields: " + getListOfUsedFields(", ") + ".");
    else
        logfile->list("Will not print any VCF " + _type + " fields.");
}


//TVCFInfoFields
TVCFInfoFields::TVCFInfoFields():TVCFFieldVector("INFO"){
    add("DP", "Number=1,Type=Integer,Description=\"Total Depth\"");
}

//TVCFGenotypeFields
TVCFGenotypeFields::TVCFGenotypeFields():TVCFFieldVector("FORMAT"){
    add("GT", "Number=1,Type=String,Description=\"Genotype\"");
    add("DP", "Number=1,Type=Integer,Description=\"Total Depth\"");
    add("GQ", "Number=1,Type=Integer,Description=\"Genotype quality\"");
    add("AD", "Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\"");
    add("AP", "Number=4,Type=Integer,Description=\"Phred-scaled allelic posterior probabilities for the four alleles A, C, G and T\"");
    add("GL", "Number=G,Type=Float,Description=\"Normalized genotype likelihoods\"");
    add("PL", "Number=G,Type=Integer,Description=\"Phred-scaled normalized genotype likelihoods\"");
    add("GP", "Number=G,Type=Integer,Description=\"Genotype posterior probabilities (phred-scaled)\"");
    add("AB", "Number=1,Type=Float,Description=\"Allelic imbalance\"");
    add("AI", "Number=1,Type=Float,Description=\"Binomial probability of allelic imbalance if Hz site\"");
}
