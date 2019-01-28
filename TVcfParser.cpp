	/*
 * TVcfParser.cpp
 *
 *  Created on: Jun 15, 2011
 *      Author: wegmannd
 */
#include "TVcfParser.h"
//--------------------------------------------------------------------
TVcfHeaderLine::TVcfHeaderLine(std::string & Line){
	init();
	std::string line=Line;
	int pp=line.find_first_of('<');
	line.erase(0, pp+1);
	line=extractBefore(line, '>');
	while(!line.empty()){
		std::string temp=extractBefore(line, ',');
		if(temp.find("=\"")>0){
			//contains quotes. Remove quotes and check for an even occurrence
			int numq=0;
			pp=temp.find_first_of('"');
			while(pp>=0){
				++numq;
				temp.erase(pp,1);
				pp=temp.find_first_of('"');
			}
			if(numq==1 || numq==3) temp=temp+extractBefore(line, '"');
			line.erase(0,1);
		}

		//get tag
		std::string tag=extractBefore(temp, '=');
		temp.erase(0,1);

		//check tag
		if(tag=="ID") id=temp;
		else if(tag=="Number"){
			numberString=temp;
			if(temp==".") number=99999;
			//else number=temp.toInt();
			else number=88888;
		}
		else if(tag=="Type"){
			typeString=temp;
			type=getTypefromString(typeString);
			if(type==UNKNOWN) "Error when parsing vcf header, unknown 'Type' in line '"+ Line +"'!";
		}
		else if(tag=="Description") desc=temp;
		else if(tag=="##INFO") continue; //HACK to fix an error in a previous version of Atlas.
		else throw "Error when parsing vcf header, unknown tag '" + tag + "' in line '"+ Line +"'!";
	}

	//throw error
	if(type!=FLAG && number<1) "Error when parsing vcf header, unknown 'Number' in line '"+ Line +"'!";
	if(id.empty() || number==-1 || type==UNKNOWN || desc.empty()){
		std::string errorMessage="Error when parsing vcf header, missing tag in line '" + Line + "':";
		if(id.empty()) errorMessage+=" id is empty!";
		if(number < 1) errorMessage+=" number is not a number!";
		if(type==UNKNOWN) errorMessage+=" unknown type!";
		if(desc.empty()) errorMessage+=" description is empty!";
		throw errorMessage;
	}

}

TVcfHeaderLine::TVcfHeaderLine(std::string & ID, std::string & Number, VCF_TYPE & Type, std::string & Desc){
	id=ID;
	numberString=Number;
	if(numberString==".") number=99999;
	else number=stringToInt(Number);
	type=Type;
	typeString=getStringfromType(type);
	desc=Desc;
	if(type!=FLAG && number<1) "Error when creating new vcf header line, unknown 'Number' entry '"+ Number +"'!";
}

void TVcfHeaderLine::update(std::string & Number, VCF_TYPE & Type, std::string & Desc){
	numberString=Number;
	if(numberString==".") number=99999;
	else number=stringToInt(Number);
	type=Type;
	typeString=getStringfromType(type);
	desc=Desc;
	if(type!=FLAG && number<1) "Error when creating new vcf header line, unknown 'Number' entry '"+ Number +"'!";
}

VCF_TYPE TVcfHeaderLine::getTypefromString(std::string & s){
	if(s=="Integer") return INTEGER;
	else if(s=="Float") return FLOAT;
	else if(s=="Flag") return FLAG;
	else if(s=="Character") return CHAR;
	else if(s=="String") return STRING;
	else return UNKNOWN;
}

std::string TVcfHeaderLine::getStringfromType(VCF_TYPE & type){
	if(type==INTEGER) return "Integer";
	else if(type==FLOAT) return "Float";
	else if(type==FLAG) return "Flag";
	else if(type==CHAR) return "Character";
	else if(type==STRING) return "String";
	else return "?";
}

void TVcfHeaderLine::init(){
	id="";
	number=-1;
	type=UNKNOWN;
	desc="";
}
std::string TVcfHeaderLine::getString(){
	return "<ID=" + id + ",Number=" + numberString + ",Type=" + typeString + ",Description=\"" + desc + "\">";
}
//--------------------------------------------------------------------
/*
TVcfFilter::TVcfFilter(std::string filter, long* CurrentLine){
	currentLine=CurrentLine;
	//has the format tag<val or tag>val or tag.sub<val or tag.sub>val
	if(filter.contains('>')) larger=true;
	else if(filter.contains('<')) larger=false;
	else throw "Filter '" + filter + "' is missing a '>' or '<' sign!";

	//get value
	std::string temp;
	if(larger) temp=filter.extract_after('>');
	else temp=filter.extract_after('<');
	if(!temp.isNumber()) throw "In filter '" + filter + "', the value '"+temp+"' is not a number!";
	val=temp.toDouble();

	//get field
	if(larger) temp=filter.extract_before('>');
	else temp=filter.extract_before('<');
	if(temp.contains('.')){
		sub=true;
		tag=temp.extract_before('.');
		subTag=temp.extract_after('.');
	} else {
		sub=false;
		tag=temp;
	}
};

bool TVcfFilter::pass(TVcfFormat* format, std::vector<std::string>* data){
	std::string buf=data->at(format->getCol(tag));
	float d;
	if(sub){
		//parse buf. Format is sub=data,sub2=data2 ...
		std::string temp;
		bool found=false;
		while(!buf.empty()){
			temp=buf.extract_sub_str(',');
			if(temp.contains('=') && temp.extract_before('=')==subTag){
				found=true;
				break;
			}
		}
		if(!found) throw "Sub-tag '"+subTag+"' missing in VCF file on line " +(std::string) *currentLine + "!";
		d=temp.extract_after('=').toDouble();
	} else d=buf.toDouble();

	if(larger && d>val) return true;
	else if (!larger && d<val) return true;
	return false;
}

void TVcfFilter::print(){
	if(sub) cout << "    -> Skipping calls where entry " << subTag << " in column " << tag;
	else cout << "    -> Skipping calls where column " << tag;
	if(larger) cout << " > ";
	else cout << " < ";
	cout << val << endl;
}
*/
//--------------------------------------------------------------------
/*
void TVcfSample::filter(TVcfFilter* filter){
	if(!missing){
		if(!filter->pass(format, &data))
			missing=false;
	}
}
*/
//--------------------------------------------------------------------
TVcfLine::TVcfLine(){
	positionParsed=false;  variantParsed=false;  idParsed=false;  filterParsed=false;  qualityParsed=false;  infoParsed=false; formatParsed=false; samplesParsed=false;
	lineNumber = -1;
	pos = -1;
}

TVcfLine::TVcfLine(std::string & line, unsigned int & numCols, long & LineNumber){
	update(line, numCols, LineNumber);
}

void TVcfLine::update(std::string & line, unsigned int & numCols, long & LineNumber){
	positionParsed=false;  variantParsed=false;  idParsed=false;  filterParsed=false;  qualityParsed=false;  infoParsed=false; formatParsed=false; samplesParsed=false;
	lineNumber=LineNumber;
	pos=-1;
	chr.clear();
	variants.clear();
	info.clear();
	formatOrdered.clear();
	format.clear();
	samples.clear();
	id.clear();
	qual.clear();
	filter.clear();

	//now read new data
	trimString(line);
	fillVectorFromStringWhiteSpace(line, data);
	if(data.size() != numCols) throw "Wrong number of columns (" + toString(data.size()) + " instead of " + toString(numCols) + ") in VCF file on line " + toString(lineNumber) + "!";
}

bool TVcfLine::variantExists(char & var){
	for(std::vector<char>::iterator it=variants.begin(); it!=variants.end(); ++it){
		if(*it==var) return true;
	}
	return false;
}

bool TVcfLine::addVariant(char var){
	if(variantExists(var)) return false;
	else variants.push_back(var);
	return true;
}

void TVcfLine::writeVariant(std::ostream & out){
	std::vector<char>::iterator it=variants.begin();
	out << *it << "\t";
	if(variants.size()>1){
		++it; out << *it; ++it;
		while(it!=variants.end()){
			out << "," << *it;
			++it;
		}
	} else out << ".";
}
//--------------------------------------------------------------------
void TVcfParser::parsePosition(TVcfLine & line){
	if(!line.positionParsed){
		//remove 'chr' and turn into int, if possible
	  /*int x=line.data[cols->Chr].find("chr");
		if(x==0) line.data[cols->Chr].remove(x, 3);
		int chr_int=line.data[cols->Chr].toInt();
		if(chr_int>0) line.data[cols->Chr]=chr_int; */
		//if(line.chr<=0) throw "Unknown chromosome '" + line.data[cols->Chr] + "' in VCF file on line " +toString(line.lineNumber) + "!";

		//just use string
		line.chr=line.data[cols->Chr];
		line.pos=stringToLong(line.data[cols->Pos]);
		if(line.pos<=0) throw "Unknown position '" + line.data[cols->Pos] + "' in VCF file on line " +toString(line.lineNumber) + "!";
		line.positionParsed=true;
	}
}

void TVcfParser::parseVariant(TVcfLine & line){
	if(!line.variantParsed){
//		std::cout << "parsing " << line.data[cols->Pos] << std::endl;
		//parse reference bases
		if(line.data[cols->Ref].length() ==1){
			line.variants.push_back(line.data[cols->Ref].c_str()[0]);
			if(line.variants[0]!='A' && line.variants[0]!='G' && line.variants[0]!='C' && line.variants[0]!='T' && line.variants[0]!='N') throw "Unknown reference allele '" + line.data[cols->Ref] + "' in VCF file on line " + toString(line.lineNumber) + "!";
		}else{
			line.variants.push_back('D'); // for deletions
			//throw "Unknown reference allele '" + line.data[cols->Ref] + "' in VCF file on line " +toString(line.lineNumber) + "!";
		}

		//alternative bases can be a comma separated list
		std::string buf;
		char var;
		if(line.data[cols->Alt]!="."){ //only if there are alternative bases
			while(!line.data[cols->Alt].empty()){
				buf=extractBefore(line.data[cols->Alt], ',');
				line.data[cols->Alt].erase(0,1);
				if(buf.length() == 1){
					var=buf.c_str()[0];
					if(var!='A' && var!='G' && var!='C' && var!='T') throw "Unknown alternative allele '" + toString(var) + "' in VCF file on line " + toString(line.lineNumber) + "!";
				}else{
					if(buf=="<NON_REF>") var='X';
					else var= 'I'; //insertion
				}
				if(!line.addVariant(var)){
					if(var != 'I' && var != 'X'){
	                    throw (std::string) "Allele '" + var + "' given multiple times in VCF file on line " + toString(line.lineNumber) + "!";
	                }
	            }
			}
		}
		line.variantParsed=true;
	}


	/*			if(buf=="<NON_REF>"){
					var='X';
				} else{
					if(buf.length() !=1){
						var='I'; // for insertions
						//throw "Unknown alternative allele '" + buf + "' in VCF file on line " + toString(line.lineNumber) + "!";
					}else{
						var=buf.c_str()[0];
						if(var!='A' && var!='G' && var!='C' && var!='T') throw "Unknown alternative allele '" + toString(var) + "' in VCF file on line " + toString(line.lineNumber) + "!";
					}
				}
				if(!line.addVariant(var)) throw "Allele '" + toString(var) + "' given multiple times in VCF file on line " + toString(line.lineNumber) + "!";
			}
		}
		line.variantParsed=true;
	}*/
}

void TVcfParser::addInfo(std::string & Line){
	TVcfHeaderLine l(Line);
	info[l.id]=l;
}
void TVcfParser::updateInfo(std::string ID, std::string Number, VCF_TYPE Type, std::string Desc){
	//check if id exists
	std::map<std::string, TVcfHeaderLine>::iterator it=info.find(ID);
	if(it==info.end()){
		//add new info
		TVcfHeaderLine l(ID, Number, Type, Desc);
		info[l.id]=l;
	} else {
		//update
		it->second.update(Number, Type, Desc);
	}
}
void TVcfParser::addFormat(std::string & Line){
	TVcfHeaderLine l(Line);
	format[l.id]=l;
}

void TVcfParser::addSample(std::string & Name){
	samples.push_back(Name);
	maxIndColPlusOne=cols->FirstInd+samples.size();
}

void TVcfParser::updateInfo(TVcfLine & line, std::string & Id, std::string & Data){
	//update or add?
	std::map<std::string, std::vector<std::string> >::iterator it=line.info.find(Id);
	if(it==line.info.end()){
		//add new tag
		if(info.find(Id)==info.end()) throw "Can not modify info, unknown info entry '" + Id + "'!";
		line.info[Id]=std::vector<std::string>();
		it=line.info.find(Id);
		it->second.push_back(Data);
	} else {
		//replace what is there
		it->second.clear();
		it->second.push_back(Data);
	}
}

void TVcfParser::updatePL(TVcfLine & line, std::string & Data, unsigned int & sample){
	//update or add?
	int col=getFormatCol(line, "PL");
	if(col < 0) throw "Did not find 'PL' info field!";
	line.samples[sample].data[col] = Data;
}

void TVcfParser::addToInfo(TVcfLine & line, std::string & Id, std::string & Data){
	//update or add?
	std::map<std::string, std::vector<std::string> >::iterator it=line.info.find(Id);
	if(it==line.info.end()){
		//add new tag
		if(info.find(Id)==info.end()) throw "Can not modify info, unknown info entry '" + Id + "'!";
		line.info[Id]=std::vector<std::string>();
		it=line.info.find(Id);
		it->second.push_back(Data);
	} else {
		//add data, if not yet there
		bool exists;
		for(std::vector<std::string>::iterator i=it->second.begin(); i!=it->second.end(); ++i){
			if((*i)==Data){
				exists=true;
				break;
			}
		}
		if(!exists) it->second.push_back(Data);
	}
}

GTLikelihoods TVcfParser::genotypeLikelihoods(TVcfLine & line, unsigned int & s){
	if(s >= line.samples.size()) throw "Sample " + toString(s) + " does not exists!";
	GTLikelihoods gt;
	if(line.samples[s].missing){
		gt.AA=-1; gt.AB=-1; gt.BB=-1;
	} else {
		int col=getFormatCol(line, "PL");
		std::string d=line.samples[s].data[col];
		//std::string a=d.extract_sub_str(',');
		stringToDouble(extractBefore(d, ','));
		gt.AA=dePhred(stringToDouble(extractBefore(d, ',')));
		d.erase(0,1);
		gt.AB=dePhred(stringToDouble(extractBefore(d, ',')));
		d.erase(0,1);
		gt.BB=dePhred(stringToDouble(d));
	}
	return gt;
}

GTLikelihoods TVcfParser::genotypeLikelihoodsPhred(TVcfLine & line, unsigned int & s){
	if(s >= line.samples.size()) throw "Sample " + toString(s) + " does not exists!";
	GTLikelihoods gt;
	if(line.samples[s].missing){
		gt.AA=-1; gt.AB=-1; gt.BB=-1;
	} else {
		int col=getFormatCol(line, "PL");
		std::string d=line.samples[s].data[col];
		gt.AA=stringToDouble(extractBefore(d, ','));
		d.erase(0,1);
		gt.AB=stringToDouble(extractBefore(d, ','));
		d.erase(0,1);
		gt.BB=stringToDouble(d);
	}
	return gt;
}

void TVcfParser::fillGenotypeLikelihoods(TVcfLine & line, unsigned int & s, double* gtl){
	if(s >= line.samples.size()) throw "Sample " + toString(s) + " does not exists!";
	if(line.samples[s].missing){
		gtl[0] = 1.0; gtl[1] = 1.0; gtl[2] = 1.0;
	} else {
		int col=getFormatCol(line, "PL");
		std::string d = line.samples[s].data[col];
		//std::string a=d.extract_sub_str(',');
		//stringToDouble(extractBefore(d, ','));
		gtl[0] = dePhred(stringToDouble(extractBefore(d, ',')));
		d.erase(0,1);
		gtl[1] = dePhred(stringToDouble(extractBefore(d, ',')));
		d.erase(0,1);
		gtl[2] = dePhred(stringToDouble(extractBefore(d, ',')));
	}
}

std::string TVcfParser::sampleContentAt(TVcfLine & line, std::string & tag, unsigned int & sample){
	checkSampleNum(line, sample);
	int col=getFormatCol(tag, line);
	if(col<0) throw "Column '"+tag+"' is missing at position " + toString(line.pos) + " on " + line.chr + "!";
	return line.samples[sample].data[col];
}

std::string TVcfParser::sampleContentAtNoCheckForMissingSample(TVcfLine & line, std::string & tag, unsigned int & sample){
//	checkSampleNum(line, sample);
	int col=getFormatCol(tag, line);
	if(col<0) return "0";//throw "Column '"+tag+"' is missing at position " + toString(line.pos) + " on " + line.chr + "!";
	return line.samples[sample].data[col];
}

int TVcfParser::getSampleNum(std::string & Name){
	int col=0;
	for(std::vector<std::string>::iterator it=samples.begin(); it!=samples.end(); ++it, ++col){
		if(*it == Name) return col;
	}
	throw "Sample '" + Name + "' is missing in the vcf file!";
	return -1;
}

std::string TVcfParser::getSampleName(unsigned int & sample){
	if(sample < 0 || sample > samples.size()) throw "Sample " + toString(sample) + " does not exists!";
	return samples[sample];
}

int TVcfParser::getNumSamples(){
	return samples.size();
}

//--------------------------------------------------------------------
//get variant info
//--------------------------------------------------------------------
std::string TVcfParser::getChr(TVcfLine & line){
	if(!line.positionParsed){
		//cerr << endl << endl << endl << "THIS LINE:" << endl << line.lineNumber << endl;
		throw "Position has not been parsed!";
	}
	return line.chr;
}

long TVcfParser::getPos(TVcfLine & line){
	if(!line.positionParsed){
		//cerr << endl << endl << endl << "THIS LINE:" << endl << line.lineNumber << endl;
		throw "Position has not been parsed!";
	}
	return line.pos;
}

int TVcfParser::getNumAlleles(TVcfLine & line){
	if(!line.variantParsed){
		throw "Position has not been parsed!";
	}
	return line.variants.size();
}

char TVcfParser::getRefAllele(TVcfLine & line){
	if(!line.variantParsed){
		throw "Position has not been parsed!";
	}
	return line.variants[0];
}

char TVcfParser::getFirstAltAllele(TVcfLine & line){
	if(!line.variantParsed){
		throw "Position has not been parsed!";
	}
	return line.variants[1];
}

char TVcfParser::getAllele(TVcfLine & line, int num){
	if(!line.variantParsed){
		throw "Position has not been parsed!";
	}
	return line.variants[num];
}

//--------------------------------------------------------------------

void TVcfParser::checkSampleNum(TVcfLine & line, unsigned int & sample){
	if(sample >= line.samples.size()) throw "Sample " + toString(sample) + " does not exists!";
	if(line.samples[sample].missing) throw "Sample " + toString(sample) + " is missing!";
}

void TVcfParser::addInfoToSample(TVcfLine & line, unsigned int & sample, std::string & tag, std::string & Data){
	if(sample >= line.samples.size()) throw "Sample " + toString(sample) + " does not exists!";
	if(!line.samples[sample].missing){
		//find position in format string
		int col=addFormatCol(tag, line);
		line.samples[sample].updateData(col, Data);
	}
}

void TVcfParser::setSampleMissing(TVcfLine & line, unsigned int & sample){
	if(sample >= line.samples.size()) throw "Sample " + toString(sample) + " does not exists!";
	line.samples[sample].missing=true;
}
void TVcfParser::setSampleHasUndefinedGenotype(TVcfLine & line, unsigned int & sample){
	if(sample >= line.samples.size()) throw "Sample " + toString(sample) + " does not exists!";
	line.samples[sample].unknownGenotype=true;
}

void TVcfParser::updateField(TVcfLine & line, std::string & tag, std::string & Data, unsigned int & sample){
	checkSampleNum(line, sample);
	int col = getFormatCol(tag, line);
	if(col < 0) throw "Column '" + tag + "' is missing at position " + toString(line.pos) + " on " + line.chr + "!";
	line.samples[sample].data[col] = Data;
};

bool TVcfParser::sampleIsHomoRef(TVcfLine & line, unsigned int & sample){
	checkSampleNum(line, sample);
	if(line.samples[sample].genotype.first==0 && line.samples.at(sample).genotype.second==0) return true;
	return false;
}

bool TVcfParser::sampleIsHeteroRefNonref(TVcfLine & line, unsigned int & sample){
	checkSampleNum(line, sample);
	if(line.samples[sample].genotype.first==0 && line.samples.at(sample).genotype.second!=0) return true;
	if(line.samples[sample].genotype.first!=0 && line.samples.at(sample).genotype.second==0) return true;
	return false;
}

char TVcfParser::getFirstAlleleOfSample(TVcfLine & line, const unsigned int & sample){
	return line.variants[line.samples[sample].genotype.first];
}

char TVcfParser::getSecondAlleleOfSample(TVcfLine & line, const unsigned int & sample){
	return line.variants[line.samples[sample].genotype.second];
}

short TVcfParser::sampleGenotype(TVcfLine & line, const unsigned int & sample){
	//NOTE: only makes sense for biallelic sites! Missing = 3
	if(line.samples[sample].missing)
		return 3;
	return line.samples[sample].genotype.first + line.samples[sample].genotype.second;
}

bool TVcfParser::sampleIsMissing(TVcfLine & line, unsigned int & s){
	if(s >= line.samples.size()) throw "Sample " + toString(s) + " does not exists!";
	return line.samples[s].missing;
}

bool TVcfParser::sampleHasUndefinedGenotype(TVcfLine & line, unsigned int & s){
	if(s >= line.samples.size()) throw "Sample " + toString(s) + " does not exists!";
	return line.samples[s].unknownGenotype;
}

float TVcfParser::sampleGenotypeQuality(TVcfLine & line, unsigned int & sample){
//	std::cerr << "check quality of sample " << sample << " and position " << line.pos << ": " << std::flush;
	checkSampleNum(line, sample);
	int col=getFormatCol(line, "GQ");
//	std::cerr << " col=" << col << std::flush;
	if(col<0) col=getFormatCol(line, "RGQ");
	if(col<0) throw "Column 'GQ' is missing at position " + toString(line.pos) + " on " + line.chr + "!";
//	std::cerr << " qual=" << line.samples[sample].data[col] << std::endl;
	return stringToDouble(line.samples[sample].data[col]);
}

int TVcfParser::phred(double x){
	return (double) -10.0 * log10(x);
}
double TVcfParser::dePhred(double x){
	return pow(10.0, -x/10.0);
}

void TVcfParser::parseFormat(TVcfLine & line){
	if(!line.formatParsed){
		std::string buf;
		int i=0;
		while(!line.data[cols->Format].empty()){
			buf=extractBefore(line.data[cols->Format], ':');
			trimString(buf);
			line.data[cols->Format].erase(0,1);

			line.format.insert(std::pair<std::string, int>(buf, i));
			line.formatOrdered.push_back(buf);
			++i;
		}
		line.formatParsed=true;
	}
}

int TVcfParser::getFormatCol(std::string & tag, TVcfLine & line){
	std::map<std::string, int>::iterator it=line.format.find(tag);
	if(it==line.format.end()) return -1;
	return it->second;
}

bool TVcfParser::formatColExists(std::string & tag, TVcfLine & line){
	std::map<std::string, int>::iterator it=line.format.find(tag);
	if(it == line.format.end()) return false;
	return true;
}

int TVcfParser::addFormatCol(std::string & tag, TVcfLine & line){
	//col exists?
	int col=getFormatCol(tag, line);
	if(col<0){
		//col does not exists -> add!
		col=line.format.size();
		line.format.insert(std::pair<std::string, int>(tag, col));
		line.formatOrdered.push_back(tag);
		//add emtpy string to all samples
		for(lineSampleIt=line.samples.begin(); lineSampleIt!=line.samples.end(); ++lineSampleIt){
			if(!(*lineSampleIt).missing) (*lineSampleIt).addData("");
		}
	}
	return col;
}

void TVcfParser::parseInfo(TVcfLine & line){
	if(!line.infoParsed){
		std::string buf, temp;
		std::map<std::string, std::vector<std::string> >::iterator it;
		while(!line.data[cols->Info].empty()){
			buf=extractBefore(line.data[cols->Info], ';');
			trimString(buf);
			line.data[cols->Info].erase(0,1);

			temp=extractBefore(buf, '=');
			buf.erase(0,1);

			line.info[temp]=std::vector<std::string>();
			it=line.info.find(temp);
			while(!buf.empty()){
				it->second.push_back(extractBefore(buf, ','));
				buf.erase(0,1);
			}
		}
		line.infoParsed=true;
	}
}

void TVcfParser::parseSamples(TVcfLine & line){
	if(!line.samplesParsed){

		parseVariant(line);
		parseFormat(line);

		//do it for all samples
		int gtCol = getFormatCol(genotypeTag, line);
		std::string gt;
		for(int i=cols->FirstInd; i<maxIndColPlusOne; ++i){
			line.samples.push_back(TVcfSample());
			//parse into std::vector (split by ':')
			while(!line.data[i].empty() > 0){
				line.samples.rbegin()->addData(extractBefore(line.data[i], ':'));
				line.data[i].erase(0,1);
			}
			//get genotype
			if(gtCol >= 0){
				gt = line.samples.rbegin()->getCol(gtCol);
				//check if data is missing: GT is either "." or "./." or ".|."
				if(gt=="." || gt=="./." || gt==".|."){
					line.samples.rbegin()->missing=true;
				} else {
					line.samples.rbegin()->missing=false;
					line.samples.rbegin()->hasGenotype=true;
					//check genotype
					if(gt.length()!=3 || (gt[1]!='/' && gt[1]!='|')) throw  "Unknown genotype '"+gt+"' in VCF file on line " + toString(line.lineNumber) + "!\n";
				}
				if(gt[1]=='/'){
					line.samples.rbegin()->setGenotype(stringToInt(extractBefore(gt, '/')), stringToInt(extractAfter(gt, '/')));
				}else if(gt[1]=='|'){
					line.samples.rbegin()->setGenotype(stringToInt(extractBefore(gt, '|')), stringToInt(extractAfter(gt, '|')));
				}
				if(!line.samples.rbegin()->checkGenotype(line.variants.size()-1)){
					std::cerr << "Unknown genotype '"+gt+"' in VCF file on line " + toString(line.lineNumber) + "!\n";
					line.samples.rbegin()->unknownGenotype = true;
				}
			} else {
				line.samples.rbegin()->missing=false;
				line.samples.rbegin()->hasGenotype=false;
			}
		}
		line.samplesParsed=true;
	}
}

void TVcfParser::writeColumnDescriptionHeader(std::ostream & out){
	out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	//samples
	for(std::vector<std::string>::iterator it=samples.begin(); it!=samples.end(); ++it){
		out << "\t" << *it;
	}
	out << std::endl;
}
void TVcfParser::writeInfoHeader(std::ostream & out){
	for(std::map<std::string, TVcfHeaderLine>::iterator it=info.begin(); it!=info.end(); ++it){
		out << "##INFO=" << it->second.getString() << std::endl;
	}
}
void TVcfParser::writeFormatHeader(std::ostream & out){
	for(std::map<std::string, TVcfHeaderLine>::iterator it=format.begin(); it!=format.end(); ++it){
		out << "##FORMAT=" << it->second.getString() << std::endl;
	}
}

void TVcfParser::writeLine(TVcfLine & line, std::ostream & out){
	//position
	if(line.positionParsed) out << line.chr << "\t" << line.pos;
	else out << line.data[cols->Chr] << "\t" << line.data[cols->Pos];

	//id
	if(line.idParsed) out << "\tERROR"; // we do not yet parse id
	else out << "\t" << line.data[cols->Id];

	//variant
	if(line.variantParsed){
		out << "\t";
		line.writeVariant(out);
	} else out << "\t" << line.data[cols->Ref] << "\t" << line.data[cols->Alt];

	//qual
	if(line.qualityParsed) out << "\tERROR"; // we do not yet parse id
	else out << "\t" << line.data[cols->Qual];

	//filter
	if(line.filterParsed) out << "\tERROR"; // we do not yet parse id
	else out << "\t" << line.data[cols->Filter];

	//info
	if(line.infoParsed){
		out << "\t";
		std::vector<std::string>::iterator i;
		bool first=true;
		for(std::map<std::string, std::vector<std::string> >::iterator it=line.info.begin(); it!=line.info.end(); ++it){
			if(first) first=false;
			else out << ";";
			out << it->first;
			if(it->second.size()>0){
				out << "=" << it->second[0];
				for(i=it->second.begin()+1; i!=it->second.end(); ++i){
					out << "," << *i;
				}
			}
		}
	} else out << "\t" << line.data[cols->Info];

	//format
	if(line.formatParsed){
		std::vector<std::string>::iterator it=line.formatOrdered.begin();
		out << "\t" << *it; ++it;
		for(;it!=line.formatOrdered.end(); ++it) out << ":" << *it;
	} else out << "\t" << line.data[cols->Format];

	//samples
	if(line.samplesParsed){
		for(std::vector<TVcfSample>::iterator it=line.samples.begin(); it!=line.samples.end(); ++it){
			it->write(out, line.formatOrdered.size());
		}
	} else {
		for(int i=cols->FirstInd; i<maxIndColPlusOne; ++i){
			out << "\t" << line.data[i];
		}
	}

	out << std::endl;
}







