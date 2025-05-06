#include "TSimulatorReference.h"
#include "coretools/Main/TLog.h"
#include "coretools/Strings/stringConstants.h"
#include "coretools/Strings/stringManipulations.h"

namespace Simulations {
using coretools::instances::logfile;

TSimulatorReference::TSimulatorReference(std::string_view Filename) : _writer(Filename, 70) {
	// open FASTA file for reference sequences
	logfile().list("Will write reference sequence to '", Filename, "'.");
}

TSimulatorReference::~TSimulatorReference() {
	_writeRefToFasta();
}

void TSimulatorReference::_writeRefToFasta() {
	if (_ref.empty()) return;
	_writer.write(_ref);
	_ref.clear();
}

void TSimulatorReference::setChr(std::string_view ChrName, long ChrLength) {
	// write if not yet written
	_writeRefToFasta();
	_writer.newContig(ChrName);

	// move to new chr
	_ref.resize(ChrLength);
}
}
