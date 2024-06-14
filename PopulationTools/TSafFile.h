/*
 * TSafFile.h
 *
 *  Created on: 2023-06-26
 *      Author: Andreas
 */

#ifndef TSAFFILE_H_
#define TSAFFILE_H_

#include "coretools/Files/TStdWriter.h"
#include "coretools/traits.h"

#include "TBgzWriter.h"

namespace PopulationTools {

class TSafFile {
	std::string _prefix;
	std::string _chr;

	uint64_t _offsetFreq = 0;
	uint64_t _offsetPos  = 0;
	size_t _nSites       = 0;
	size_t _sumBand      = 0;

	GLF::TBGzWriter _freqWriter;
	GLF::TBGzWriter _posWriter;
	coretools::TStdWriter _idxWriter;

	void _writeIdx() {
		_idxWriter.write(_chr.size());
		_idxWriter.write(_chr);
		_idxWriter.write(_nSites);
		_idxWriter.write(_sumBand);

		_idxWriter.write(_offsetPos);
		_idxWriter.write(_offsetFreq);
	}

public:
	TSafFile(std::string_view Prefix, size_t nSamples)
		: _prefix(Prefix), _freqWriter(_prefix + ".saf.gz"),
		  _posWriter(_prefix + ".saf.pos.gz"), _idxWriter(_prefix + ".saf.idx") {
		constexpr char magic[8] = "safv4";
		_freqWriter.write(magic);
		_posWriter.write(magic);

		_idxWriter.write(magic);
		_idxWriter.write(static_cast<size_t>(2 * nSamples)); // assume diploid

		_offsetFreq = _freqWriter.tell();
		_offsetPos  = _posWriter.tell();
	}


	~TSafFile() {
		if (_nSites > 0) _writeIdx();
	}
	TSafFile(const TSafFile &)            = delete;
	TSafFile &operator=(const TSafFile &) = delete;
	TSafFile(TSafFile &&)                 = default;
	TSafFile &operator=(TSafFile &&)      = default;

	template<typename LogContainer>
	void write(std::string_view Chr, size_t Pos, const LogContainer& AlleleFreqs, size_t lower=0) {
		if (_chr.empty()) { // first chromosome
			_chr = Chr;
		} else if (_chr != Chr) { // change of chromosome
			_writeIdx();

			// reset
			_chr        = Chr;
			_nSites     = 0;
			_sumBand    = 0;
			_offsetFreq = _freqWriter.tell();
			_offsetPos  = _posWriter.tell();
		}

		static std::vector<float> floats;
		floats.clear();

		for (size_t i = lower; i <AlleleFreqs.size(); ++i) {
			floats.push_back(static_cast<float>(coretools::underlying(AlleleFreqs[i])));
		}

		_freqWriter.write(static_cast<int>(lower));
		_freqWriter.write(static_cast<int>(floats.size()));
		_freqWriter.write(floats);

		_posWriter.write(static_cast<int>(Pos)); // 0-indexed in saf-file

		++_nSites;
		_sumBand += floats.size();
	}
};

} // namespace PopulationTools

#endif
