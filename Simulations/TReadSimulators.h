/*
 * TReadSimulator.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TREADSIMULATORS_H_
#define SIMULATIONS_TREADSIMULATORS_H_

#include "TReadSimulator.h"
#include "../GenotypeLikelihoods/SequencingError/TModels.h"
#include "TPostMortemDamage.h"
#include "TGenomePosition.h"

namespace Simulations {


class TReadSimulators {
private:
	double _averageReadLength = 0;
	double _maxFragmentLength = 0;

	// simulation tools
	BAM::TReadGroups _readGroups;
	GenotypeLikelihoods::TPostMortemDamage _PMD;
	GenotypeLikelihoods::SequencingError::TModels _recal;

	// read simulator
	std::vector<std::unique_ptr<TReadSimulator>> _readSimulators;
	std::vector<coretools::Probability> _simGroupFrequencies;
	std::vector<coretools::Probability> _cumulSimGroupFrequenies;

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
	void writeUnwrittenAlignments(const genometools::TGenomePosition & Position, TSimulatorBamFile &BamFile);
	void simulate(const genometools::TGenomePosition & Position, const std::vector<Base>& Haplotype, TSimulatorBamFile &BamFile);

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
