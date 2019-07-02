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

    TVCFField(const std::string & Tag, const std::string & Description);

    ~TVCFField(){}

    void writeHeader(std::string type, gz::ogzstream & vcf);
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
    TVCFFieldVector(std::string Type);

    std::string type();
    void add(std::string tag, std::string description);
    TVCFField* getFieldPointer(const std::string & tag);
    bool fieldExists(const std::string & tag);
    void acceptField(const std::string & tag);
    bool useField(std::string tag);
    void clearUsed();
    int numUsed();
    std::string getListOfUsedFields(std::string delim);
    void fillVectorWithTagsOfUsedFields(std::vector<std::string> & vec);
    void writeVCFHeader(gz::ogzstream & vcf);
	//report VCF fields
    void reportUsedFields(TLog* logfile);
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
// ADD NEW FIELDS BELOW!
////////////////////////////////////////////////////////////////////////////////////////////////////////

class TVCFInfoFields:public TVCFFieldVector{
public:
    TVCFInfoFields();
};

class TVCFGenotypeFields:public TVCFFieldVector{
public:
    TVCFGenotypeFields();;
};




#endif /* TVCFFIELDS_H_ */
