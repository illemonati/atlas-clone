/*
 * TGlfMultiReader.h
 *
 *  Created on: Feb 11, 2020
 *      Author: wegmannd
 */

#ifndef GLF_TGLFMULTIREADER_H_
#define GLF_TGLFMULTIREADER_H_

#include <string>
#include <vector>

#include "TGLF.h"

#include "genometools/GenotypeTypes.h"
#include "genometools/TFastaReader.h"

using TGenotypeLikelihoodsOneAllelicCombinationVector = std::vector<genometools::TGenotypeLikelihoodsOneAllelicCombination>;
using TGenotypeLikelihoodsAllCombinationsVector                      = std::vector<genometools::TGenotypeLikelihoodsAllCombinations>;

namespace GLF {

//----------------------------------------------------
// TGlfMultiReader
//----------------------------------------------------
class TGlfMultiReader {
private:
	size_t _numGLFs = 0;
	std::vector<std::string> _GLFNames;
	std::vector<TGlfReader> _GLFs;
	bool _readersOpened = false;

	// active files
	// Object will loop only over active files
	size_t _minSamplesWithData  = 0;
	std::vector<bool> _GLFIsActive;
	std::vector<TGlfReader *> _activeGLFs;

	// Moving along active files
	TGlfChromosome _curChr;
	uint32_t _curRefId = 0;
	uint32_t _minDepth = 0;
	size_t _windowStart = 0;
	size_t _windowSize  = 64;
	std::vector<TGenotypeLikelihoodsAllCombinationsVector> _dataWindow;
	std::vector<size_t> _numActive;

	// reference
	genometools::TFastaReader fastaReader;

	void _openGLFs();
	int _getGLFIndexFromName(const std::string &name) const;
	void _setActive(size_t index);
	void _setAllInactive();
	void _prepareParsing();
	bool _jumpToNextPosition();

	bool _moveToNextChromosome();

public:
	const TGenotypeLikelihoodsAllCombinationsVector& data(size_t iWindow) const noexcept {return _dataWindow[iWindow];};

	TGlfMultiReader();

	void openGLFs(const std::vector<std::string> &Filenames);
	void openGLFs();
	void addReference(const std::string& FastaFile);
	const genometools::TChromosomes& chromosomes();	
	void minSamplesWithData(size_t MinSamplesWithData) { _minSamplesWithData = MinSamplesWithData; };

	// set active / inactive
	void setActive(int index);
	void setActive(const std::string &name);
	void setActive(int index1, int index2);
	void setActive(const std::string &name1, const std::string &name2);
	void setActive(const std::vector<int> &indexes);
	void setActive(const std::vector<std::string> &names);
	void setAllActive();

	// parse
	std::vector<size_t> readWindow();

	// output
	std::vector<std::string> namesOfActiveFiles() const;
	std::vector<std::string> sampleNamesOfActiveFiles() const;

	// access data
	uint32_t numActiveSamples() const noexcept { return _activeGLFs.size(); }
	std::string chr() const { return _curChr.name(); }
	constexpr uint32_t position(size_t iWindow) const noexcept { return _windowStart + iWindow; }
	bool hasRef() const noexcept {return fastaReader.isOpen();}
	coretools::TView<genometools::Base> refView() const {
		assert(hasRef());
		return fastaReader.view(_curRefId, _windowStart, std::min(_windowSize, _curChr.length() - _windowStart));
	}
};

}; // end namespace GLF

#endif /* GLF_TGLFMULTIREADER_H_ */
