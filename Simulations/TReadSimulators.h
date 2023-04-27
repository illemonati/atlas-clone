/*
 * TReadSimulator.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TREADSIMULATORS_H_
#define SIMULATIONS_TREADSIMULATORS_H_

#include "TReadSimulator.h"
#include "TPostMortemDamage.h"
#include "genometools/GenomePositions/TGenomePosition.h"
#include "../BAM/TBamFile.h"

namespace Simulations {


class TReadSimulators {
private:
	double _averageReadLength = 0;
	double _maxFragmentLength = 0;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;

	// read simulator
	std::vector<std::unique_ptr<TReadSimulator>> _readSimulators;
	std::vector<coretools::Probability> _simGroupFrequencies;
	std::vector<coretools::Probability> _cumulSimGroupFrequencies;

	// function to initialize read groups
	void _initializeReadGroups(const BAM::RGInfo::TReadGroupInfo & RGinfo);
	void _initializeContamination(bool &perReadGroup, std::map<std::string, double> &contaminationMap);
	void _initializeReadSimulator();
	void _initializeReadGroupFrequencies(const BAM::RGInfo::TReadGroupInfo & RGinfo);
	void _determineMaxFragmentLength();

public:
	TReadSimulators(const std::string & RgInfoFileName);
	TReadSimulators(TReadSimulators && other) = default;

	//interact
	void simulate(const genometools::TGenomePosition & Position, const std::vector<Base>& Haplotype, Simulations::TSimulatedOutputFile &SimFile);

	//getters
	[[nodiscard]] std::unique_ptr<TReadSimulator>& sample();
	[[nodiscard]] double maxFragmentLength() const { return _maxFragmentLength; };
	[[nodiscard]] double averageFragmentLength() const { return _averageReadLength; };
	BAM::TReadGroups& readGroups() { return _readGroups; };

	std::vector<std::unique_ptr<TReadSimulator>>::iterator begin(){ return _readSimulators.begin(); };
	std::vector<std::unique_ptr<TReadSimulator>>::iterator end(){ return _readSimulators.end(); };
};

} // end namespace Simulations

#endif /* SIMULATIONS_TREADSIMULATORS_H_ */
