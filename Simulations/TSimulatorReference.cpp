#include "TSimulatorReference.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringConstants.h"
#include "coretools/Strings/stringManipulations.h"

namespace Simulations {
using coretools::instances::logfile;

TSimulatorReference::TSimulatorReference(std::string_view Filename): _filename(Filename) {
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
	_writeRefToFasta();
	_fasta.close();
	_fastaIndex.close();
}

void TSimulatorReference::_writeRefToFasta() {
	if (_ref.empty()) return;
	// write to fasta
	_fasta << ">" << _chrName;
	for (size_t l = 0; l < _ref.size(); ++l) {
		if (l % 70 == 0) _fasta << "\n";
		_fasta << _ref[l];
	}
	_fasta << "\n";

	// add to index
	_oldOffset += _chrName.size() + 2;
	_fastaIndex << coretools::str::readBefore(_chrName, coretools::str::whitespaces) << "\t" << _ref.size() << "\t" << _oldOffset
				<< "\t70\t71\n";
	_oldOffset += _ref.size() + (int)(_ref.size() / 70);
	if (_ref.size() % 70 != 0) _oldOffset += 1;
	_ref.clear();
}

void TSimulatorReference::setChr(std::string_view ChrName, long ChrLength) {
	// write if not yet written
	_writeRefToFasta();

	// move to new chr
	_chrName   = ChrName;
	_ref.resize(ChrLength);
}
}
