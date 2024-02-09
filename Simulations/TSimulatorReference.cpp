#include "TSimulatorReference.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringManipulations.h"

namespace Simulations {
using coretools::instances::logfile;

TSimulatorReference::TSimulatorReference(std::string Filename){
	open(Filename);
}

void TSimulatorReference::open(std::string Filename){
	_filename = Filename;
	// open FASTA file for reference sequences
	logfile().list("Will write reference sequence to '" + _filename + "'.");
	_fasta.open(_filename.c_str());
	if (!_fasta) UERROR("Failed to open file '", _filename, "' for writing!");
	_filename += ".fai";
	_fastaIndex.open(_filename.c_str());
	if (!_fastaIndex) UERROR("Failed to open file '", _filename, "' for writing!");
	_oldOffset = 0;
}

TSimulatorReference::~TSimulatorReference() {
	if (_chrName != "" && _needsWriting) _writeRefToFasta();
	_fasta.close();
	_fastaIndex.close();
}

void TSimulatorReference::_writeRefToFasta() {
	// write to fasta
	_fasta << ">" << _chrName;
	for (int l = 0; l < _chrLength; ++l) {
		if (l % 70 == 0) _fasta << "\n";
		_fasta << _ref[l];
	}
	_fasta << "\n";

	// add to index
	std::string tmp = _chrName;
	_oldOffset += _chrName.size() + 2;
	_fastaIndex << coretools::str::extractBeforeWhiteSpace(tmp) << "\t" << _chrLength << "\t" << _oldOffset
				<< "\t70\t71\n";
	_oldOffset += _chrLength + (int)(_chrLength / 70);
	if (_chrLength % 70 != 0) _oldOffset += 1;

	_needsWriting = false;
}

void TSimulatorReference::_allocateStorage(long length) {
	// allocate storage
	_ref.resize(length);
}

void TSimulatorReference::setChr(std::string ChrName, long ChrLength) {
	if(!_fasta){
		DEVERROR("Fasta file not opened yet!");
	}
	// write if not yet written
	if (_chrName != "" && _needsWriting) _writeRefToFasta();

	// move to new chr
	_chrName   = ChrName;
	_chrLength = ChrLength;
	_ref.resize(_chrLength);
	_needsWriting = true;
}
}
