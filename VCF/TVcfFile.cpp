/*
 * TVcfFile.cpp
 *
 *  Created on: Aug 8, 2011
 *      Author: wegmannd
 */

#include "TVcfFile.h"

namespace VCF{

//--------------------------------------------------------------------
TVcfFile_base::TVcfFile_base(std::string & Filename, bool zipped){
	filename = Filename;
	eof=false;
	openStream(filename, zipped);
};

TVcfFile_base::TVcfFile_base(){
	currentLine=0;
	automaticallyWriteVcf=false;
	eof=false;
	numCols=-1;
	totalFileSize=-1;
	myOutStream=NULL;
	myStream=NULL;
	inputStreamOpend=false;
	outputStreamOpend=false;
};

void TVcfFile_base::openStream(const std::string & filename){
	if(coretools::str::readAfterLast(filename, '.') == "gz"){
		openStream(filename, true);
	} else {
		openStream(filename, false);
	}
};

void TVcfFile_base::openStream(const std::string & Filename, const bool & zipped){
	//open stream
	filename = Filename;
	if(zipped) myStream = new gz::igzstream(filename.c_str());
	else myStream = new std::ifstream(filename.c_str());
	if(!(*myStream) || myStream->fail() || !myStream->good()) throw "Failed to open file '" + filename + "'!";
	inputStreamOpend = true;

	//first line contains fileformat version
	std::string temp;
	getline(*myStream, temp);
	currentLine=1;

	if(stringContains(temp, "##fileformat")){
		fileFormat = coretools::str::extractAfter(temp, '=');
		coretools::str::trimString(fileFormat);
		if(fileFormat!="VCFv4.0" && fileFormat!="VCFv4.1" && fileFormat!="VCFv4.2" && fileFormat!="VCFv4.3") throw "VCF file is not in 'VCFv4.0' format!";
	} else throw "Missing VCF file format specification on first line! Is '" + Filename + "' a VCf file?";

	//parse header
	parseHeaderVCF_4_0();
};

void TVcfFile_base::openOutputStream(const std::string & Filename, const bool & zipped){
	//open stream
	outputFilename = Filename;
	if(zipped) myOutStream = new gz::ogzstream(outputFilename.c_str());
	else myOutStream = new std::ofstream(outputFilename.c_str());
	if(!(*myOutStream)) throw "Failed to open file '" + outputFilename + "'!";
	outputStreamOpend = true;
};

void TVcfFile_base::setOutStream(std::ostream & os){
	myOutStream = &os;
}

void TVcfFile_base::parseHeaderVCF_4_0(){
	//read the header of the vcf file and stop after that
	std::string temp, buf;
	bool headerRowRead=false;
	while(!myStream->eof() && myStream->peek()=='#'){
		++currentLine;
		getline(*myStream, temp);
		if(stringContains(temp, "#CHROM")){
			if(headerRowRead) throw "Found more than one header row!";
			//analyze header: save which column contains the chromosome, position, refbase, altbases, info, format and species
			coretools::str::trimString(temp);
			int i=0;
			while(!temp.empty()){
				buf = coretools::str::extractBeforeWhiteSpace(temp);
				coretools::str::trimString(buf);
				temp.erase(0,1);
				if(i<parser.cols.FirstInd) parser.cols.set(buf, i);
				else parser.addSample(buf);
				++i;
			}
			numCols=i;
			headerRowRead=true;
		} else {
			//parse all other header lines
			if(temp.find("##FORMAT")==0) parser.addFormat(temp);
			else if(temp.find("##INFO")==0) parser.addInfo(temp);
			else unknownHeader.push_back(temp);
		}
	}
	//check cols
	parser.cols.check();
	if(parser.samples.size()<1) throw "VCF file contains no samples!";
	currentLine = 0;
}

void TVcfFile_base::addNewHeaderLine(std::string headerLine){
	unknownHeader.push_back(headerLine);
}

/*
void TVcfFile::addFilter(my_string filter){
	filters.push_back(TVcfFilter(filter, &currentLine));
	applyFilters=true;
}

void TVcfFile::filterSamples(){
	if(!samplesParsed) throw TException("Unable to filter samples: samples not parsed yet!");
	if(applyFilters){
		for(vector<TVcfFilter>::iterator itF=filters.begin(); itF!=filters.end(); ++itF){
			for(vector<TVcfSample>::iterator it=samples.begin(); it!=samples.end(); ++it)
				it->filter(&(*itF));
		}
	}
}

void TVcfFile::printFilters(){
	if(applyFilters){
		cout << " - Using the following filters:" << endl;
		for(vector<TVcfFilter>::iterator itF=filters.begin(); itF!=filters.end(); ++itF)
			itF->print();
	}
}
*/

GTLikelihoods TVcfFile_base::_genotypeLikelihoods(TVcfLine* line, unsigned int s){
	return parser.genotypeLikelihoods(*line, s);
}

GTLikelihoods TVcfFile_base::_genotypeLikelihoodsPhred(TVcfLine* line, unsigned int s){
	return parser.genotypeLikelihoodsPhred(*line, s);
}

void TVcfFile_base::fillGenotypeLiklihoods(TVcfLine* line, unsigned int sample, float* gtl){
	parser.fillGenotypeLikelihoods(*line, sample, gtl);
}

void TVcfFile_base::fillPherdScore(TVcfLine* line, unsigned int sample, uint8_t & gtl_0, uint8_t & gtl_1, uint8_t & gtl_2){
	parser.fillPhredScore(*line, sample, gtl_0, gtl_1, gtl_2);
}

void TVcfFile_base::fillLog10GenotypeLikelihoods(TVcfLine* line, unsigned int sample, double & gtl_0, double & gtl_1, double & gtl_2){
	parser.fillLog10GenotypeLikelihoods(*line, sample, gtl_0, gtl_1, gtl_2);
}

int TVcfFile_base::sampleNumber(std::string & Name){
	return parser.getSampleNum(Name);
}

std::string TVcfFile_base::sampleName(unsigned int num){
	return parser.getSampleName(num);
}

int TVcfFile_base::numSamples(){
	return parser.getNumSamples();
}

bool TVcfFile_base::sampleIsMissing(TVcfLine* line, unsigned int sample){
	return parser.sampleIsMissing(*line, sample);
}

bool TVcfFile_base::sampleHasUnknownGenotype(TVcfLine* line, unsigned int sample){
	return parser.sampleHasUndefinedGenotype(*line, sample);
}

std::string TVcfFile_base::fieldContentAsString(std::string tag, TVcfLine* line, unsigned int sample){
	return parser.sampleContentAt(*line, tag, sample);
}

int TVcfFile_base::fieldContentAsInt(std::string tag, TVcfLine* line, unsigned int sample){
	return coretools::str::convertString<int>(parser.sampleContentAt(*line, tag, sample));
}

int TVcfFile_base::depthAsIntNoCheckForMissingSample(std::string tag, TVcfLine* line, unsigned int sample){
	return coretools::str::convertString<int>(parser.sampleContentAtNoCheckForMissingSample(*line, tag, sample));
}

void TVcfFile_base::setSampleMissing(TVcfLine* line, unsigned int sample){
	parser.setSampleMissing(*line, sample);
}

void TVcfFile_base::setSampleHasUndefinedGenotype(TVcfLine* line, unsigned int sample){
	parser.setSampleHasUndefinedGenotype(*line, sample);
}

void TVcfFile_base::updateField(TVcfLine* line, std::string & tag, std::string & Data, unsigned int sample){
	parser.updateField(*line, tag, Data, sample);
};

void TVcfFile_base::writeHeaderVCF_4_0(){
	*myOutStream << "##fileformat=" << fileFormat << std::endl;
	parser.writeFormatHeader(*myOutStream);
	parser.writeInfoHeader(*myOutStream);
	//write unknown header columns
	for(std::vector<std::string>::iterator i=unknownHeader.begin(); i!=unknownHeader.end(); ++i){
		*myOutStream << *i << std::endl;
	}
	//add column description
	parser.writeColumnDescriptionHeader(*myOutStream);
}

bool TVcfFile_base::readLine(){
	std::string temp;
	do{
		if(myStream->eof()){
			eof=true;
			return false;
		}
		++currentLine;
		getline(*myStream, temp);
	} while(temp.empty());
	tempLine.update(temp, numCols, currentLine);

	//what to parse?
	for(std::vector<pt2Function>::iterator it=usedParsers.begin(); it!=usedParsers.end(); ++it){
		(parser.*(*it))(tempLine);
	}
	return true;
}

void TVcfFile_base::updateInfo(TVcfLine* line, std::string id, std::string data){

}
void TVcfFile_base::addToInfo(TVcfLine* line, std::string Id, std::string Data){

}

//--------------------------------------------------------------------
TVcfFileSingleLine::TVcfFileSingleLine(std::string & filename, bool zipped){
	openStream(filename, zipped);
	written = true;
	eof = false;
}

TVcfFileSingleLine::~TVcfFileSingleLine(){
	if(!written){
		if(automaticallyWriteVcf) parser.writeLine(tempLine, *myOutStream);
	}

}

void TVcfFileSingleLine::writeLine(){
	parser.writeLine(tempLine, *myOutStream);
	written = true;
}

bool TVcfFileSingleLine::next(){
	//write old Line
	if(!written){
		if(automaticallyWriteVcf) writeLine();
	}
	if(!readLine()) return false;
	else {
		written=false;
		return true;
	}
}

void TVcfFileSingleLine::updateInfo(std::string Id, std::string Data){
	if(written) throw "Can not update line, no line read!";
	parser.updateInfo(tempLine, Id, Data);
}
void TVcfFileSingleLine::updatePL(std::string Data, unsigned int sample){
	parser.updatePL(tempLine, Data, sample);
}
void TVcfFileSingleLine::addToInfo(std::string Id, std::string Data){
	if(written) throw "Can not update line, no line read!";
	parser.addToInfo(tempLine, Id, Data);
}
std::string TVcfFileSingleLine::fieldContentAsString(std::string tag, unsigned int sample){
	return TVcfFile_base::fieldContentAsString(tag, &tempLine, sample);
}
int TVcfFileSingleLine::fieldContentAsInt(std::string tag, unsigned int sample){
	return TVcfFile_base::fieldContentAsInt(tag, &tempLine, sample);
}

int TVcfFileSingleLine::depthAsIntNoCheckForMissingSample(std::string tag, unsigned int sample){
	return TVcfFile_base::depthAsIntNoCheckForMissingSample(tag, &tempLine, sample);
}

GTLikelihoods TVcfFileSingleLine::genotypeLikelihoods(unsigned int sample){
	return _genotypeLikelihoods(&tempLine ,sample);
}
GTLikelihoods TVcfFileSingleLine::genotypeLikelihoodsPhred(unsigned int sample){
	return _genotypeLikelihoodsPhred(&tempLine ,sample);
}
void TVcfFileSingleLine::fillGenotypeLikelihoods(unsigned int sample, float* gtl){
	fillGenotypeLiklihoods(&tempLine, sample, gtl);
}

void TVcfFileSingleLine::fillPhredScore(unsigned int sample, uint8_t & gtl_0, uint8_t & gtl_1, uint8_t & gtl_2){
	fillPherdScore(&tempLine, sample, gtl_0, gtl_1, gtl_2);
}

void TVcfFileSingleLine::fillLog10GenotypeLikelihoods(unsigned int sample, double & gtl_0, double & gtl_1, double & gtl_2){
	TVcfFile_base::fillLog10GenotypeLikelihoods(&tempLine, sample, gtl_0, gtl_1, gtl_2);
}

uint64_t TVcfFileSingleLine::position(){
	return parser.getPos(tempLine);
}
uint64_t TVcfFileSingleLine::positionZeroBased(){
	return position() - 1;
}
std::string TVcfFileSingleLine::chr(){
	return parser.getChr(tempLine);
}
int TVcfFileSingleLine::getNumAlleles(){
	return parser.getNumAlleles(tempLine);
}
std::string TVcfFileSingleLine::getRefAllele(){
	return parser.getRefAllele(tempLine);
}
std::string TVcfFileSingleLine::getFirstAltAllele(){
	return parser.getFirstAltAllele(tempLine);
}
std::string TVcfFileSingleLine::getAllele(int num){
	return parser.getAllele(tempLine, num);
};

bool TVcfFileSingleLine::isBialleleicSNP(){
	if(parser.getNumAlleles(tempLine)==2){
		std::string ref = parser.getAllele(tempLine, 0);
		std::string alt = parser.getAllele(tempLine, 1);
		return (ref=="A" || ref=="C" || ref=="G" || ref=="T") && (alt=="A" || alt=="C" || alt=="G" || alt=="T");
	}
	return false;
};

bool TVcfFileSingleLine::variantQualityIsMissing(){
	return parser.variantQualityIsMissing(tempLine);
}

double TVcfFileSingleLine::variantQuality(){
	return parser.variantQuality(tempLine);
}

void TVcfFileSingleLine::setSampleMissing(unsigned int sample){
	TVcfFile_base::setSampleMissing(&tempLine, sample);
}
void TVcfFileSingleLine::setSampleHasUndefinedGenotype(unsigned int sample){
	TVcfFile_base::setSampleHasUndefinedGenotype(&tempLine, sample);
}
void TVcfFileSingleLine::updateField(std::string tag, std::string & Data, unsigned int sample){
	TVcfFile_base::updateField(&tempLine, tag, Data, sample);
}
bool TVcfFileSingleLine::sampleIsMissing(unsigned int sample){
	return TVcfFile_base::sampleIsMissing(&tempLine, sample);
}
bool TVcfFileSingleLine::sampleHasUndefinedGenotype(unsigned int sample){
	return TVcfFile_base::sampleIsMissing(&tempLine, sample);
}
bool TVcfFileSingleLine::sampleIsHaploid(unsigned int sample){
	return parser.sampleIsHaploid(tempLine, sample);
};
bool TVcfFileSingleLine::sampleIsDiploid(unsigned int sample){
	return parser.sampleIsDiploid(tempLine, sample);
};
bool TVcfFileSingleLine::sampleIsHomoRef(unsigned int sample){
	return parser.sampleIsHomoRef(tempLine, sample);
}
bool TVcfFileSingleLine::sampleIsHeteroRefNonref(unsigned int sample){
	return parser.sampleIsHeteroRefNonref(tempLine, sample);
}
float TVcfFileSingleLine::sampleGenotypeQuality(unsigned int sample){
	return parser.sampleGenotypeQuality(tempLine, sample);
}
BAM::Base TVcfFileSingleLine::getFirstAlleleOfSample(unsigned int num){
	return parser.getFirstAlleleOfSample(tempLine, num)[0];
};

BAM::Base TVcfFileSingleLine::getSecondAlleleOfSample(unsigned int num){
	return parser.getSecondAlleleOfSample(tempLine, num)[0];
};
BAM::Genotype TVcfFileSingleLine::sampleGenotype(const unsigned int & num){
	return BAM::Genotype(parser.getFirstAlleleOfSample(tempLine, num)[0], parser.getSecondAlleleOfSample(tempLine, num)[0]);
};

// int TVcfFileSingleLine::sampleDepth(unsigned int sample){
double TVcfFileSingleLine::sampleDepth(unsigned int sample){
	//return 0 if sample is missing
	if(parser.sampleIsMissing(tempLine, sample))
		return 0;
	//check if depth is given
	std::string DP = "DP";
	if(parser.formatColExists(DP, tempLine))
		return coretools::str::convertString<double>(parser.sampleContentAt(tempLine, DP, sample));
		// return convertString<int>(parser.sampleContentAt(tempLine, DP, sample));
	else return -1;
}

}; //end namespace

