/*
 * TReadSimulator.h
 *
 *  Created on: Aug 25, 2022
 *      Author: phaentu
 */

#ifndef SIMULATIONS_TREADSIMULATORS_H_
#define SIMULATIONS_TREADSIMULATORS_H_

#include "TReadSimulator.h"

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
	void _initializePMD();
	void _initializeQualityTransformations();
	void _initializeContamination(bool &perReadGroup, std::map<std::string, double> &contaminationMap);
	void _initializeReadSimulator();
	void _initializeReadGroupFrequencies(const BAM::RGInfo::TReadGroupInfo & RGinfo);
	void _determineMaxFragmentLength();

public:
	TReadSimulators(const std::string & RgInfoFileName);

	//getters
	[[nodiscard]] TReadSimulator& sample();
	[[nodiscard]] double maxFragmentLength() const { return _maxFragmentLength; };
	[[nodiscard]] double averageFragmentLength() const { return _averageReadLength; };


};

} // end namespace Simulations

#endif /* SIMULATIONS_TREADSIMULATORS_H_ */
