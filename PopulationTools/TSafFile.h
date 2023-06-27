/*
 * TSafFile.h
 *
 *  Created on: 2023-06-26
 *      Author: Andreas
 */

#ifndef TSAFFILE_H_
#define TSAFFILE_H_

#include "coretools/Containers/TView.h"
#include "coretools/Files/TWriterImpl.h"
#include "coretools/devtools.h"
#include "coretools/traits.h"
#include <type_traits>
namespace PopulationTools {

class TSafFile {
	std::string _prefix;
	std::string _chr;

	uint64_t _offsetFreq = 0;
	uint64_t _offsetPos  = 0;
	size_t _nSites       = 0;
	size_t _sumBand      = 0;

	coretools::TBGzWriter _freqWriter;
	coretools::TBGzWriter _posWriter;
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
	TSafFile(std::string_view Prefix, size_t nSamples) :_prefix(Prefix), _freqWriter(_prefix + ".saf.gz"), _posWriter(_prefix + ".saf.pos.gz"), _idxWriter(_prefix + ".saf.idx") {
		constexpr char magic[8] = "safv4";
		_freqWriter.write(magic);
		_posWriter.write(magic);

		_idxWriter.write(magic);
		_idxWriter.write(static_cast<size_t>(2*nSamples)); // assume diploid

		_offsetFreq = _freqWriter.tell();
		_offsetPos = _posWriter.tell();
	}

	~TSafFile() {
		if (_nSites > 0) _writeIdx();
	}
	TSafFile(const TSafFile &)            = delete;
	TSafFile &operator=(const TSafFile &) = delete;
	TSafFile(TSafFile &&)                 = default;
	TSafFile &operator=(TSafFile &&)      = default;

	template<typename T>
	void write(std::string_view Chr, size_t Pos, const std::vector<T> &AlleleFreqs) {
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

		_freqWriter.write(static_cast<int>(0)); // We're not doing this
		_freqWriter.write(static_cast<int>(AlleleFreqs.size()));
		if constexpr (std::is_same_v<T, float>) {
			_freqWriter.write(AlleleFreqs);
		} else {
			std::vector<float> floats;
			floats.reserve(AlleleFreqs.size());
			for (auto af : AlleleFreqs) { floats.push_back(static_cast<float>(coretools::underlying(af))); }
			_freqWriter.write(floats);
		}

		_posWriter.write(static_cast<int>(Pos - 1)); // 0-indexed in saf-file

		++_nSites;
		_sumBand += AlleleFreqs.size();
	}
};

} // namespace PopulationTools

#endif
