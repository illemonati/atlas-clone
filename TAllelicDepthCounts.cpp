/*
 * TAllelicDepthCounts.cpp
 *
 *  Created on: Feb 10, 2020
 *      Author: wegmannd
 */

#include "TAllelicDepthCounts.h"

void TAllelicDepthCounts::_allocateStorage(const uint32_t MaxAllelicDepth){
	_maxAllelicDepth = MaxAllelicDepth;
	_size = _maxAllelicDepth + 1;

	//allocate
	_counts = new uint32_t***[_size];
	for(uint32_t i=0; i<_size; ++i){
		_counts[i] = new uint32_t**[_size];
		for(uint32_t j=0; j<_size; ++j){
			_counts[i][j] = new uint32_t*[_size];
			for(uint32_t k=0; k<_size; ++k){
				_counts[i][j][k] = new uint32_t[_size];
			}
		}
	}

	clear();
};

void TAllelicDepthCounts::clear(){
	for(uint32_t i=0; i<_size; ++i){
		for(uint32_t j=0; j<_size; ++j){
			for(uint32_t k=0; k<_size; ++k){
				for(uint32_t l=0; l<_size; ++l){
					_counts[i][j][k][l] = 0;
				}
			}
		}
	}
};

void TAllelicDepthCounts::_freeStorage(){
	for(uint32_t i=0; i<_size; ++i){
		for(uint32_t j=0; j<_size; ++j){
			for(uint32_t k=0; k<_size; ++k){
				delete[] _counts[i][j][k];
			}
			delete[] _counts[i][j];
		}
		delete[] _counts[i];
	}
	delete[] _counts;
};


void TAllelicDepthCounts::addSite(uint32_t numA, uint32_t numC, uint32_t numG, uint32_t numT){
	if(numA < _size && numC < _size && numG < _size && numT < _size)
		++_counts[numA][numC][numG][numT];
};


void TAllelicDepthCounts::addSiteZeroDepth(){
	++_counts[0][0][0][0];
};

void TAllelicDepthCounts::write(const std::string filename, bool printEmpty){
	//open file
	TOutputFileZipped out(filename);

	//write header
	out.writeHeader({"A", "C", "G", "T", "Depth", "majorAllele", "minorAllele", "Counts"});

	//write counts
	//uint32_t max = 0;
	for(uint32_t i=0; i<_size; ++i){
		for(uint32_t j=0; j<_size; ++j){
			for(uint32_t k=0; k<_size; ++k){
				for(uint32_t l=0; l<_size; ++l){
					if(printEmpty || _counts[i][j][k][l] > 0){
						//write numA, C, G and T and depth
						out << i << j << k << l << i+j+k+l;

						//find max
						uint32_t max = i;
						if(j > max) max = j;
						if(k > max) max = k;
						if(l > max) max = l;

						//identify those that are at max
						std::vector<char> tmp;
						if(i == max) tmp.push_back('A');
						if(j == max) tmp.push_back('C');
						if(k == max) tmp.push_back('G');
						if(l == max) tmp.push_back('T');

						//write major
						out << tmp[0];

						//find minor
						if(tmp.size() > 1){
							out << tmp[1];
						} else {
							//find second
							uint32_t second = 0;
							if(i < max && i > second) second = i;
							if(j < max && j > second) second = j;
							if(k < max && k > second) second = k;
							if(l < max && l > second) second = l;

							//print minor
							if(i == second) out << 'A';
							else if(j == second) out << 'C';
							else if(k == second) out << 'G';
							else out << 'T';
						}

						//write counts
						out << _counts[i][j][k][l] << std::endl;
					}
				}
			}
		}
	}
};




