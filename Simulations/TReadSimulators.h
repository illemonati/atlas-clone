/*
 * TReadSimulator.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TREADSIMULATORS_H_
#define SIMULATIONS_TREADSIMULATORS_H_

#include "coretools/Main/TRandomPicker.h"
#include "genometools/GenomePositions/TGenomePosition.h"

#include "PMD/TModels.h"
#include "TReadSimulator.h"
#include "genometools/Genotypes/TwoBases.h"

namespace Simulations {


class TReadSimulators {
private:
	double _averageReadLength = 0;
	double _maxFragmentLength = 0;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::PMD::TModels _pmd;
	GenotypeLikelihoods::SequencingError::TModels _recal;

	// read simulator
	std::vector<std::unique_ptr<TReadSimulator>> _readSimulators;
	std::vector<coretools::Probability> _simGroupFrequencies;
	coretools::TRandomPicker _picker;

	// function to initialize read groups
	void _initializeReadGroups(const BAM::RGInfo::TReadGroupInfo & RGinfo);
	void _initializeReadSimulator();
	void _initializeReadGroupFrequencies(const BAM::RGInfo::TReadGroupInfo & RGinfo);
	void _determineMaxFragmentLength();

public:
	TReadSimulators(std::string_view RGOutName, std::string_view RGInName);
	TReadSimulators(TReadSimulators && other) = default;

	//interact
	std::pair<size_t, size_t> simulate(const genometools::TGenomePosition &Position,
									   const std::vector<genometools::TwoBase> &Haplotype,
									   coretools::TView<genometools::Base> Reference, BAM::TOutputBamFile &BamFile);

	// getters
	[[nodiscard]] double maxFragmentLength() const { return _maxFragmentLength; };
	[[nodiscard]] double averageReadLength() const { return _averageReadLength; };
	BAM::TReadGroups& readGroups() { return _readGroups; };
	const TReadSimulator& readSimulator(size_t RG) const noexcept {return *_readSimulators[RG];}
	size_t numRG() const noexcept {return _readSimulators.size();}

	std::vector<std::unique_ptr<TReadSimulator>>::iterator begin(){ return _readSimulators.begin(); };
	std::vector<std::unique_ptr<TReadSimulator>>::iterator end(){ return _readSimulators.end(); };
};

} // end namespace Simulations

#endif /* SIMULATIONS_TREADSIMULATORS_H_ */
