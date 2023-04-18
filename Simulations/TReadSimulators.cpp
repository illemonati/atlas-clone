/*
 * TReadSimulator.cpp
 *
 */

#include "TReadGroupInfo.h"
#include "coretools/Strings/stringFunctions.h"
#include "coretools/Main/TParameters.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TRandomGenerator.h"
#include "TReadSimulators.h"
#include "coretools/Types/probability.h"
#include "coretools/Main/globalConstants.h"

namespace Simulations {

using BAM::RGInfo::TReadGroupInfo;
using coretools::instances::logfile;
using coretools::instances::randomGenerator;

void TReadSimulators::_initializeReadGroups(const TReadGroupInfo & RGinfo) {
	// create simulation read groups
	using BAM::RGInfo::InfoType;
	for(size_t i = 0; i < RGinfo.size(); ++i){
		logfile().startIndent("Read group '", RGinfo[i].name(), "':");
		std::string type = RGinfo[i].getString(InfoType::seqType);
		logfile().list("Sequencing type: ", type);
		logfile().list("Frequency: ", _simGroupFrequencies[i]);

		//initialize by type
		if(type == "single"){
			_readSimulators.push_back(std::make_unique<TReadSimulatorSingleEnd>(_readGroups[i], RGinfo[i]));
		} else if(type == "paired"){
			_readSimulators.push_back(std::make_unique<TReadSimulatorPairedEnd>(_readGroups[i], RGinfo[i]));
		} else {
			UERROR("Unable to understand read group type '" + type + "'! Use either 'single' or 'paired'.");
		}
		logfile().endIndent();
	}
}

void TReadSimulators::_initializeReadGroupFrequencies(const TReadGroupInfo &RGinfo) {
	_cumulSimGroupFrequenies.resize(RGinfo.size());
	_simGroupFrequencies.resize(RGinfo.size());

	using BAM::RGInfo::InfoType;
	// fill frequencies and cumulative frequencies
	std::vector<double> tmp;
	RGinfo.fillContainerPerReadGroup(tmp, InfoType::RGFrequency);
	coretools::fillFromNormalized(_simGroupFrequencies, tmp);
	coretools::fillCumulative(_simGroupFrequencies, _cumulSimGroupFrequenies);
}

void TReadSimulators::_determineMaxFragmentLength(){
	// precalculate some stuff
	_averageReadLength = 0;
	_maxFragmentLength = 0;

	for (size_t i = 0; i < _readSimulators.size(); ++i) {
		_averageReadLength += _simGroupFrequencies[i] * _readSimulators[i]->meanReadLength();
		if (_readSimulators[i]->maxFragmentLength() > _maxFragmentLength){
			_maxFragmentLength = _readSimulators[i]->maxFragmentLength();
		}
	}

	if(_averageReadLength < 1.0){
		UERROR("Chosen parameters result in an average fragment length across read groups < 1.0!");
	}
}

TReadSimulators::TReadSimulators(const std::string & RgInfoFileName){
	// Read sequencing parameters from RG Info / Command line
	TReadGroupInfo RGinfo;
	_readGroups = RGinfo.createReadGroups(RgInfoFileName);
	using BAM::RGInfo::InfoType;
	logfile().addIndent();
	RGinfo.parse(InfoType::seqType, InfoType::cycles, InfoType::fragmentLength, InfoType::baseQuality, InfoType::mappingQuality, InfoType::softClipping, InfoType::pmd, InfoType::recal, InfoType::RGFrequency);
	logfile().endIndent();

	// complete RG details
	for (auto &rg : _readGroups) {
		rg.sequencingCenter_CN =
			coretools::__GLOBAL_APPLICATION_NAME__ + " " + coretools::__GLOBAL_APPLICATION_VERSION__;
		rg.description_DS          = "Simulated with commit " + coretools::__GLOBAL_APPLICATION_COMMIT__;
		rg.sequencingTechnology_PL = "ILLUMINA";
	}

	_initializeReadGroupFrequencies(RGinfo);

	//Initialize read groups
	logfile().startIndent("Will use the following ", _readGroups.size(), " read groups:");
	_initializeReadGroups(RGinfo);

	// B) initialize PMD
	//------------------
	if (coretools::instances::parameters().parameterExists("pmd")) {
		_PMD.initialize(coretools::instances::parameters().getParameter("pmd"), _readGroups);
	} else {
		_PMD.initialize(RGinfo);
	}

	// add PMD to simulators
	// TODO: also initialize PMD from RGInfo
	for (size_t r = 0; r < _readSimulators.size(); ++r) {
		_readSimulators[r]->setPMD(&_PMD[r]);
	}

	// C) initialize contamination
	//----------------------------
	// TODO: Think about contamination object for both estimation and simulation

	// D) other things
	//----------------
	// initialize read group frequencies frequencies

	//warn if read group info columns were not used
	logfile().endIndent();
	RGinfo.warnAboutUnusedColumnsInFile();
	logfile().endIndent();

	//prepare simulations
	_determineMaxFragmentLength();
}

void TReadSimulators::simulate(const genometools::TGenomePosition & Position, const std::vector<Base>& Haplotype, BAM::TOutputBamFile &BamFile){
	//sample which simulator to use
	size_t thisSimulator = randomGenerator().pickOne(_cumulSimGroupFrequenies);
	_readSimulators[thisSimulator]->simulate(Position, Haplotype, BamFile);
}

std::unique_ptr<TReadSimulator>& TReadSimulators::sample(){
	return _readSimulators[randomGenerator().pickOne(_cumulSimGroupFrequenies)];
}


} // end namespace Simulations


