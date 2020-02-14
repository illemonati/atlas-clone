/*
 * TBed.h
 *
 *  Created on: Mar 16, 2017
 *      Author: phaentu
 */

#ifndef TBED_H_
#define TBED_H_


class TBedChromosome{
private:
	std::string _name;
	std::map<uint64_t, uint64_t> _windows;
	std::map<uint64_t, uint64_t>::iterator _windowIt;
	bool _parsed;

public:
	TBedChromosome(std::string & Name){
		_name = Name;
		rewind();
	};

	~TBedChromosome(){
		//delete all windows
		_windows.clear();
	};

	void addWindow(uint64_t start, uint64_t end){
		if(end <= start) throw "Window  [" + toString(start) + ", " + toString(end) + "] is not valid!";
		//check if there is overlap with any other window
		for(_windowIt=_windows.begin(); _windowIt!=_windows.end(); ++_windowIt){
			if((start >= _windowIt->first && start < _windowIt->second) || (end > _windowIt->first && end <= _windowIt->second))
				throw "Error reading bed file: window [" + toString(start) + ", " + toString(end) + ") overlaps window [" + toString(_windowIt->first) + ", " + toString(_windowIt->second) + ")!";
		}
		_windows.emplace(start, end);
	};

	void print(){
		std::cout << "Chromosome '" << _name << "':" << std::endl;
		for(_windowIt=_windows.begin(); _windowIt!=_windows.end(); ++_windowIt){
			std::cout << " - [" << _windowIt->first << ", " << _windowIt->second << "]" << std::endl;
		}
	};

	size_t size(){
		return _windows.size();
	};

	void setAsParsed(){
		_parsed = true;
	};

	void rewind(){
		_parsed = false;
	};

	bool parsed(){
		return _parsed;
	};

	bool begin(){
		_windowIt = _windows.begin();
		if(_windowIt == _windows.end()) return false;
		return true;
	};

	bool next(){
		++_windowIt;
		if(_windowIt == _windows.end()) return false;
		return true;
	};

	bool reachedEnd(){
		if(_windowIt == _windows.end()) return true;
		return false;
	};

	uint64_t curStart(){
		return _windowIt->first;
	};

	uint64_t curEnd(){
		return _windowIt->second;
	};

	uint64_t curSize(){
		return _windowIt->second - _windowIt->first;
	};
};

class TBed{
private:
	std::map<std::string, TBedChromosome*> _chromosomes;
	std::map<std::string, TBedChromosome*>::iterator _chrIt;
	std::string _curChr;
	bool _allChrParsed;

	void readFile(){
		//open file
		std::istream* myStream = NULL;
		if(filename.find(".gz")) myStream = new gz::igzstream(filename.c_str());
		else myStream = new std::ifstream(filename.c_str());
		if(!*myStream) throw "Failed to open BED file '" + filename + "'!";

		//tmp variables
		long lineNum = 0;
		std::vector<std::string> vec;
		_curChr = "";

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
				if(vec[0] != _curChr){
					_chrIt = _chromosomes.find(vec[0]);
					if(_chrIt == _chromosomes.end()){
						_chromosomes.insert(std::pair<std::string, TBedChromosome*>(vec[0], new TBedChromosome(vec[0])));
						_chrIt = _chromosomes.find(vec[0]);
					}
					_curChr = vec[0];
				}

				//add positions
				_chrIt->second->addWindow(stringToLong(vec[1]), stringToLong(vec[2]));
			}
		}

		//close file
		delete myStream;
		_chrIt = _chromosomes.end();
	};

	void _setCurChrAsParsed(){
		_chrIt->second->setAsParsed();

		//check if all chromosomes have been parsed
		_allChrParsed = true;
		for(auto& it : _chromosomes){
			if(!it.second->parsed()){
				_allChrParsed = false;
				break;
			}
		}
	};

public:
	std::string filename;

	TBed(std::string Filename){
		filename = Filename;
		readFile();
		_curChr = "";
		rewind();
	};

	~TBed(){
		//delete all chromosomes
		for(auto& it : _chromosomes){
			delete it.second;
		}
		_chromosomes.clear();
	};

	void rewind(){
		for(auto& it : _chromosomes){
			it.second->rewind();
		}

		_allChrParsed = false;
		_chrIt = _chromosomes.end();
	};

	bool setChr(const std::string & chr){
		//set cur chr as parsed
		if(_chrIt != _chromosomes.end()){

			//already on correct chromosome?
			if(_chrIt->first == chr){
				return true;
			}

			//set current as parsed
			_setCurChrAsParsed();
		}

		//jump to requested chromosome
		_curChr = chr;
		_chrIt = _chromosomes.find(_curChr);
		if(_chrIt!=_chromosomes.end()){
			_chrIt->second->begin();
			return true;
		}
		return false;
	};

	std::string curChr(){
		if(_chrIt==_chromosomes.end()) return "";
		return _chrIt->first;
	};

	bool nextWindow(){
		if(_chrIt==_chromosomes.end()) return false;
		if(!_chrIt->second->next()){
			//set current as parsed
			_setCurChrAsParsed();
			return false;
		}
		return true;
	};

	bool reachedEnd(){
		return _allChrParsed;
	};

	bool reachedEndOfChr(){
		if(_chrIt==_chromosomes.end()) return true;
		return _chrIt->second->reachedEnd();
	};

	uint64_t curWindowStart(){
		if(_chrIt==_chromosomes.end()) return -1;
		return _chrIt->second->curStart();
	};

	uint64_t curWindowEnd(){
		if(_chrIt==_chromosomes.end()) return -1;
		return _chrIt->second->curEnd();
	};

	uint64_t curWindowSize(){
		if(_chrIt==_chromosomes.end()) return 0;
				return _chrIt->second->curSize();
	};

	void print(){
		std::cout << "Bed File '" << filename << "':" << std::endl;
		for(_chrIt=_chromosomes.begin(); _chrIt!=_chromosomes.end(); ++_chrIt)
			_chrIt->second->print();
	};

	uint64_t size(){
		uint64_t s = 0;
		for(auto& it : _chromosomes){
			s += it.second->size();
		}
		return s;
	};

	int getNumChromosomes(){
		return _chromosomes.size();
	};

	int getNumWindowsOnCurChr(){
		if(_chrIt==_chromosomes.end()) return 0;
		return _chrIt->second->size();
	};
};



#endif /* TBED_H_ */
