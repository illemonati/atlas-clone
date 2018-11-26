/*
 * TVCFFields.h
 *
 *  Created on: Nov 18, 2018
 *      Author: phaentu
 */

#ifndef TVCFFIELDS_H_
#define TVCFFIELDS_H_

#include <vector>
#include "gzstream.h"
#include "TLog.h"

//------------------------------------------------------
// TVCFField
//------------------------------------------------------
class TVCFField{
public:
	std::string tag;
	std::string description;
	bool accepted;
	bool used;

	TVCFField(const std::string & Tag, const std::string & Description){
		tag = Tag;
		description = Description;
		accepted = false;
		used = false;
	};

	~TVCFField(){};

	void writeHeader(std::string type, gz::ogzstream & vcf){
		vcf << "##" << type << "=<ID=" << tag << ',' << description << ">\n";
	};
};

//------------------------------------------------------
// TVCFFieldVector
//------------------------------------------------------
class TVCFFieldVector{
private:
	std::string _type;
	std::vector<TVCFField> availableFields;
	std::vector<TVCFField*> usedFields;

public:
	TVCFFieldVector(std::string Type){
		_type = Type;
	};

	std::string type(){ return _type; };

	void add(std::string tag, std::string description){
		availableFields.emplace_back(tag, description);
	};

	TVCFField* getFieldPointer(const std::string & tag){
		for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it){
			if(tag == it->tag) return &(*it);
		}
		throw "VCF " + _type + " field '" + tag + "' is unknown!";
	};

	bool fieldExists(const std::string & tag){
		for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it){
			if(it->tag == tag)
				return true;
		}
		return false;
	};

	void acceptField(const std::string & tag){
		getFieldPointer(tag)->accepted = true;
	};

	bool useField(std::string tag){
		TVCFField* pt = getFieldPointer(tag);
		if(pt->accepted){
			if(!pt->used){
				usedFields.push_back(pt);
				pt->used = true;
			}
			return true;
		}
		return false;
	};

	void clearUsed(){
		usedFields.clear();
		for(std::vector<TVCFField>::iterator it = availableFields.begin(); it != availableFields.end(); ++it)
			it->used = false;
	};

	int numUsed(){
		return usedFields.size();
	};

	std::string getListOfUsedFields(std::string delim){
		std::string out;
		for(std::vector<TVCFField*>::iterator it = usedFields.begin(); it != usedFields.end(); ++it){
			if(out.length() > 0) out += delim;
			out += (*it)->tag;
		}
		return out;
	};

	void fillVectorWithTagsOfUsedFields(std::vector<std::string> & vec){
		vec.clear();
		for(std::vector<TVCFField*>::iterator it = usedFields.begin(); it != usedFields.end(); ++it)
			vec.push_back((*it)->tag);
	};

	void writeVCFHeader(gz::ogzstream & vcf){
		//info fields
		for(std::vector<TVCFField*>::iterator it = usedFields.begin(); it != usedFields.end(); ++it)
			(*it)->writeHeader(_type, vcf);
	};

	//report VCF fields
	void reportUsedFields(TLog* logfile){
		if(numUsed() > 0)
			logfile->list("Will print these VCF " + _type + " fields: " + getListOfUsedFields(", ") + ".");
		else
			logfile->list("Will not print any VCF " + _type + " fields.");
	};
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// ADD NEW FIELDS BELOW!
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TVCFInfoFields:public TVCFFieldVector{
public:
	TVCFInfoFields():TVCFFieldVector("INFO"){
		add("DP", "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\"");
	};
};

class TVCFGenotypeFields:public TVCFFieldVector{
public:
	TVCFGenotypeFields():TVCFFieldVector("FORMAT"){
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
	};
};




#endif /* TVCFFIELDS_H_ */
