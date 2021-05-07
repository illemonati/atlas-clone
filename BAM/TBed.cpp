
#include "TBed.h"

namespace BAM{

//-------------------------------------------------------------
// TBedChromosome
//-------------------------------------------------------------

bool operator<(const std::string & Name, const TBedChromosome & Chr){
    return Name < Chr.name;
};

//-------------------------------------------------------------
// TBed_base
// base class to store a list of windows
//-------------------------------------------------------------

std::set<TBedChromosome, std::less<>>::iterator TBed_base::addChromosome(const std::string & Chr){
	//get chromosome id
	auto it = _chromosomes.find(Chr);
	if(it == _chromosomes.end()){ // new chromosome
		_chromosomeNames[_chromosomes.size()] = Chr;
        return _chromosomes.emplace(_chromosomes.size(), Chr).first;
    }
    return it;
};

void TBed_base::addChromosomes(const TChromosomes & Chromosomes){
	for(std::vector<TChromosome>::const_iterator c = Chromosomes.cbegin(); c != Chromosomes.cend(); ++c){
		auto it = _chromosomes.emplace(c->refID(), c->name);
		if(it.second){
			//insert succeeded: add name
			_chromosomeNames[c->refID()] = c->name;
		} else {
			//insert failed as element already exists: check if refID is correct!
			if(it.first->refId != c->refID()){
				throw std::runtime_error("void TBed::addChromosomes(const TChromosomes & Chromosomes): chromosome already exists, but with different ID!");
			}
		}
	}
};

bool TBed_base::hasWindowsOnChr(const uint32_t refId){
	return _chromosomeNames.find(refId) != _chromosomeNames.end();
};

bool TBed_base::hasWindowsOnChr(const std::string& Chr){
	return _chromosomes.find(Chr) != _chromosomes.end();
};

uint32_t TBed_base::numChromosomesWithWindows() const{
    throw std::runtime_error("uint32_t TBed_base::numChromosomesWithWindows() const: function not implemented for base class!");
};

uint32_t TBed_base::getRefID(const std::string& Chr){
	auto it = _chromosomes.find(Chr);
	if(it == _chromosomes.end()){
		throw std::runtime_error("uint32_t TBed_base::getRefID(const std::string Chr): chromosome '" + Chr + "' does not exist!");
	}
	return it->refId;
};

std::string TBed_base::getChromosomeName(const uint32_t refId){
	auto it = _chromosomeNames.find(refId);
	if(it == _chromosomeNames.end()){
		throw std::runtime_error("std::string TBed_base::getChromosomeName(const uint32_t refId): ref ID '" + toString(refId) + "' does not exist!");
	}
	return it->second;
};

//-------------------------------------------------------------
// TBed
//-------------------------------------------------------------
void TBed::add(TGenomeWindow Window){
	//check if chromosome exists
	if(_chromosomeNames.find(Window.refID()) == _chromosomeNames.end()){
		throw std::runtime_error("void TBed::add(TGenomeWindow Window): chromosome with refID " + toString(Window.refID()) + " does not exist!");
	}
    //insert on chromosome that already has windows
    //Note: windows are sorted by chr and start. Windows do not overlap.
    //find first window with same position or later
	auto it = _bed.lower_bound(Window);

    //merge with existing
    //can overlap with one upstream, but not more than one as otherwise the two preceding ones should overlap.
    //can overlap with many downstream (if large)
    //move one up
    if(it != _bed.begin()){
        --it;
    }

    //incorporate downstream as long as there is overlap
    while(it != _bed.end() && Window.mergeWith(*it)){
        it = _bed.erase(it);
    }

    //insert window
    _bed.insert(Window);
};

void TBed::add(const TGenomePosition & Position){
	add(TGenomeWindow(Position));
};

void TBed::add(const uint32_t RefID, const uint32_t Pos){
	//check if chromosome exists
	add(TGenomeWindow(RefID, Pos));
};

void TBed::add(const std::string& Chr, const uint32_t Pos){
	auto it = addChromosome(Chr);
	add(TGenomeWindow(it->refId, Pos));
};

void TBed::add(const std::string& Filename){
	//add windows from file
	TInputFile in(Filename, 3);

	std::vector<std::string> vec;
	while(in.read(vec)){
		//make sure chromosome exists
		auto it = addChromosome(vec[0]);

		//parse position
		uint32_t start = convertStringCheck<uint32_t>(vec[1]);
		uint32_t end = convertStringCheck<uint32_t>(vec[2]);

		add(TGenomeWindow(it->refId, start, end));
	}
};

void TBed::add(const std::string& Filename, const TChromosomes & Chromosomes){
	//add all chromosomes
	addChromosomes(Chromosomes);

	//add windows from file
	TInputFile in(Filename, 3);

	std::vector<std::string> vec;
	while(in.read(vec)){
		//parse data
		uint32_t refId = Chromosomes.refID(vec[0]);
		uint32_t start = convertStringCheck<uint32_t>(vec[1]);
		uint32_t end = convertStringCheck<uint32_t>(vec[2]);

		add(TGenomeWindow(refId, start, end));
	}
};

void TBed::write(const std::string& Filename) const{
	TOutputFile out(Filename);
	out.noHeader(3);

	for(auto& it: _bed){
		auto c = _chromosomeNames.find(it.refID());
		if(c == _chromosomeNames.end()){
			throw std::runtime_error("void TBed::write(const std::string Filename) const: chromosome with refID " + toString(it.refID()) + " does not exist!");
		}
		out << c->second << it.fromOnChr() << it.toOnChr() << std::endl;
	}
};

uint64_t TBed::size() const{
	return _bed.size();
};

uint64_t TBed::length() const{
	uint64_t l = 0;
	for(auto& it : _bed){
		l += it.size();
	}
	return l;
};

uint32_t TBed::numChromosomesWithWindows() const{
	return _numChromosomesWithWindows(_bed);
};

bool TBed::exists(const TGenomeWindow& window) const{
	//search entry according to chr and start
	auto it = _bed.find(window);
	if(it == _bed.end()){
		return false;
	}

	//now check for end
	if(it->to() == window.to()){
		return true;
	} else {
		return false;
	}
};

//-----------------------------------------------------
// TGenomeWindowList
//-----------------------------------------------------
void TGenomeWindowList::add(const TGenomeWindow& Window){
	_list.insert(Window);
};

void TGenomeWindowList::add(const std::string Filename, const TChromosomes & Chromosomes){
	//add windows from file
	TInputFile in(Filename, 3);

	std::vector<std::string> vec;
	while(in.read(vec)){
		//parse data
		uint32_t refId = Chromosomes.refID(vec[0]);
		uint32_t start = convertStringCheck<uint32_t>(vec[1]);
		uint32_t end = convertStringCheck<uint32_t>(vec[2]);

		add(TGenomeWindow(refId, start, end));
	}
};

void TGenomeWindowList::write(const std::string Filename, const TChromosomes & Chromosomes) const{
	TOutputFile out(Filename);
	out.noHeader(3);

	for(auto& it:  _list){
		out << Chromosomes.name(it.refID()) << it.fromOnChr() << it.toOnChr() << std::endl;
	}
};

uint64_t TGenomeWindowList::size() const{
	return _list.size();
};

uint64_t TGenomeWindowList::length() const{
	uint64_t l = 0;
	for(auto& it : _list){
		l += it.size();
	}
	return l;
};

uint32_t TGenomeWindowList::numChromosomesWithWindows() const{
	return _numChromosomesWithWindows(_list);
};

bool TGenomeWindowList::exists(const TGenomeWindow& window) const{
	//search entry according to chr and start
	auto it = _list.find(window);
	if(it == _list.end()){
		return false;
	}

	//now check for end
	if(it->to() == window.to()){
		return true;
	} else {
		return false;
	}
};

bool TGenomeWindowList::hasWindowsOnChr(uint32_t refId) const{
	auto it = _list.lower_bound(TGenomePosition(refId, 0));
    if(it == _list.end() || it->refID() > refId){
		return false;
	}
	return true;
};

uint32_t TGenomeWindowList::numWindowsOnChr(const uint32_t refId) const{
	auto it = _list.lower_bound(TGenomePosition(refId, 0));
	if(it == _list.end() || it->refID() > refId){
		return 0;
	} else {
		uint32_t num = 0;
		while(it != _list.end() && it->refID() == refId){
			++num;
			++it;
		}
		return num;
	}
};

/*
//---------------------------------------
// TBedChromosome
//---------------------------------------
void TBedChromosome::_fuseNeighboringWindowsIfNeeded(TBedWindowMap::iterator & it){
	TBedWindowMap::iterator before = it; --before;
	while(it != _windows.begin() && before->second >= it->first){
		//windows overlap: fuse
		before->second = std::max(before->second, it->second);

		_windows.erase(it);
		it = before;
		--before;
	}

	TBedWindowMap::iterator after = it; ++after;
	while(after != _windows.end() && after->first <= it->second){
		//windows overlap: fuse
		it->second = std::max(it->second, after->second);

		_windows.erase(after);
		after = it; ++after;
	}
};

void TBedChromosome::addWindow(uint64_t start, uint64_t end){
	if(end <= start) throw "Window  [" + toString(start) + ", " + toString(end) + "] is not valid!";
	//check if there is overlap with any other window
	for(TBedWindowMap::iterator it = _windows.begin(); it!=_windows.end(); ++it){
		if((start >= it->first && start < it->second) || (end > it->first && end <= it->second))
			throw "Error reading bed file: window [" + toString(start) + ", " + toString(end) + ") overlaps window [" + toString(_windowIt->first) + ", " + toString(_windowIt->second) + ")!";
	}
	_windows.emplace(start, end);
};

void TBedChromosome::addOrExtendWindow(const uint64_t start, const uint64_t end){
	if(end <= start) throw "Window  [" + toString(start) + ", " + toString(end) + "] is not valid!";
	//check if there is overlap with any other window

	//check if window with same start exists
	TBedWindowMap::iterator it = _windows.find(start);

	if(it != _windows.end()){
		//window has same start: extend if necessary
		if(end > it->second){
			it->second = end;
		}
	} else {
		//create new window and check if it overlaps with neighboring windows
		it = _windows.emplace(start, end).first;
		_fuseNeighboringWindowsIfNeeded(it);
	}
};

void TBedChromosome::print(const std::string & chrName){
	std::cout << "Chromosome '" << chrName << "':" << std::endl;
	for(_windowIt=_windows.begin(); _windowIt!=_windows.end(); ++_windowIt){
		std::cout << " - [" << _windowIt->first << ", " << _windowIt->second << "]" << std::endl;
	}
};

uint64_t TBedChromosome::length(){
	uint64_t l = 0;
	for(auto& it : _windows){
		l += it.second - it.first;
	}
	return l;
};

bool TBedChromosome::begin(){
	_windowIt = _windows.begin();
	if(_windowIt == _windows.end()) return false;
	return true;
};

bool TBedChromosome::next(){
	++_windowIt;
	if(_windowIt == _windows.end()) return false;
	return true;
};

bool TBedChromosome::reachedEnd(){
	if(_windowIt == _windows.end()) return true;
	return false;
};

void TBedChromosome::write(TOutputFile & out, const std::string & chrName){
	for(_windowIt=_windows.begin(); _windowIt!=_windows.end(); ++_windowIt){
		out << chrName << _windowIt->first << _windowIt->second << std::endl;
	}
};

//---------------------------------------
// TBed
//---------------------------------------
void TBed::readFile(const std::string filename){
	//open file
	std::istream* myStream = nullptr;
	if(filename.find(".gz")){
		myStream = new gz::igzstream(filename.c_str());
	} else {
		std::ifstream* tmp = new std::ifstream(filename.c_str());
		if(!tmp->is_open())
			throw "Failed to open BED file '" + filename + "'!";
		myStream = tmp;
	}
	if(!myStream->good() || myStream->eof()) throw "Failed to open BED file '" + filename + "'!";

	//tmp variables
	long lineNum = 0;
	std::vector<std::string> vec;
	_curChr = "";

	//read file
	while(myStream->good() && !myStream->eof()){
		++lineNum;
		std::string line;
		std::getline(*myStream, line);
		fillContainerFromStringWhiteSpaceSkipEmpty(line, vec);
		//skip empty lines
		if(vec.size() > 0){
			if(vec.size() < 3) throw "Less than three columns in bed file '" + filename + "' on line " + toString(lineNum) + "!";

			//get chromosome
			if(vec[0] != _curChr){
				_chrIt = _chromosomes.find(vec[0]);
				if(_chrIt == _chromosomes.end()){
					_chromosomes.emplace(vec[0], TBedChromosome{});
					_chrIt = _chromosomes.find(vec[0]);
				}
				_curChr = vec[0];
			}

			//add positions
			_chrIt->second.addWindow(stringToLong(vec[1]), stringToLong(vec[2]));
		}
	}

	//close file
	delete myStream;
	_chrIt = _chromosomes.end();

	//check if file contains windows
	if(size() == 0){
		throw "BED file '" + filename + "' is empty (or does not exists)!";
	}
};

void TBed::_setCurChrAsParsed(){
	_chrIt->second.setAsParsed();

	//check if all chromosomes have been parsed
	_allChrParsed = true;
	for(auto& it : _chromosomes){
		if(!it.second.parsed()){
			_allChrParsed = false;
			break;
		}
	}
};

void TBed::addWindowCurChr(const uint64_t start, const uint64_t end){
	_chrIt->second.addOrExtendWindow(start, end);
};

void TBed::addSite(const uint64_t pos){
	_chrIt->second.addOrExtendWindow(pos, pos+1);
};

void TBed::rewind(){
	for(auto& it : _chromosomes){
		it.second.rewind();
	}

	_allChrParsed = false;
	_chrIt = _chromosomes.end();
};

void TBed::clear(){
	//delete all chromosomes
	_chromosomes.clear();
};

bool TBed::setChr(const std::string & chr){
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
		_chrIt->second.begin();
		return true;
	}
	return false;
};

void TBed::setChrCreateIfNew(const std::string & chr){
	if(!setChr(chr)){
		//create chromosome
		_chromosomes.emplace(chr, TBedChromosome{});
		setChr(chr);
	}
};

std::string TBed::curChr(){
	if(_chrIt==_chromosomes.end()) return "";
	return _chrIt->first;
};

bool TBed::nextWindow(){
	if(_chrIt==_chromosomes.end()) return false;
	if(!_chrIt->second.next()){
		//set current as parsed
		_setCurChrAsParsed();
		return false;
	}
	return true;
};

bool TBed::reachedEndOfChr(){
	if(_chrIt==_chromosomes.end()) return true;
	return _chrIt->second.reachedEnd();
};

uint64_t TBed::curWindowStart(){
	if(_chrIt==_chromosomes.end()) return -1;
	return _chrIt->second.curStart();
};

uint64_t TBed::curWindowEnd(){
	if(_chrIt==_chromosomes.end()) return -1;
	return _chrIt->second.curEnd();
};

uint64_t TBed::curWindowSize(){
	if(_chrIt==_chromosomes.end()) return 0;
			return _chrIt->second.curLength();
};

void TBed::print(){
	std::cout << "Bed File:" << std::endl;
	for(auto& it:  _chromosomes){
		it.second.print(it.first);
	}
};

size_t TBed::size(){
	size_t s = 0;
	for(auto& it : _chromosomes){
		s += it.second.size();
	}
	return s;
};

uint64_t TBed::length(){
	uint64_t l = 0;
	for(auto& it : _chromosomes){
		l += it.second.length();
	}
	return l;
};

int TBed::getNumWindowsOnCurChr(){
	if(_chrIt==_chromosomes.end()) return 0;
	return _chrIt->second.size();
};

void TBed::write(const std::string filename){
	TOutputFile out(filename);
	out.noHeader(3);

	for(auto& it:  _chromosomes){
		it.second.write(out, it.first);
	}
};

bool TBed::test(){
	//clear
	clear();

	//add windows to check functionalities
	//create chromosome with three windows of length 20 each
	setChrCreateIfNew("chr1");
	addWindowCurChr(100, 110);
	addWindowCurChr(105, 120);

	addWindowCurChr(210, 220);
	addWindowCurChr(200, 210);

	addWindowCurChr(310, 315);
	addWindowCurChr(300, 320);

	//create chromosome with one windows of length 220
	setChrCreateIfNew("chr2");
	addWindowCurChr(100, 110);
	addWindowCurChr(200, 210);
	addWindowCurChr(300, 310);

	addWindowCurChr(110, 200);
	addWindowCurChr(200, 320);

	//create chromosome with one windows of length 201
	setChrCreateIfNew("chr3");
	for(uint64_t i = 100; i<200; ++i){
		addSite(i);
	}
	for(uint64_t i = 300; i>199; --i){
		addSite(i);
	}

	//add an extra window to chr1
	setChrCreateIfNew("chr1");
	addWindowCurChr(500, 520);

	//now test
	if(getNumChromosomes() != 3){
		print();
		throw "BED Test failed: wrong number of chromosomes!";
	}
	if(size() != 6){
		print();
		throw "BED Test failed: wrong size!";
	}
	if(length() != 4*20 + 220 + 201){
		print();
		throw "BED Test failed: wrong length!";
	}

	write("test.bed");

	return true;
};
*/

}; // end namespace
