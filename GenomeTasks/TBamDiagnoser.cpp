/*
 * TBamDiagnoser.cpp
 *
 *  Created on: May 30, 2020
 *      Author: phaentu
 */

#include "TBamDiagnoser.h"
#include "coretools/Main/TParameters.h"

namespace GenomeTasks {

using coretools::TCountDistributionVector;
using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {
void writeHistogram(const std::vector<TCountDistributionVector<>> &distVec, std::string_view outName,
					std::string_view header, const std::vector<std::string> &ReadGroupNames, TCountDistributionVector<>* all = nullptr) {
	// displays distributions of type 'TCountDistributionVector' as a histogram
	const bool hasAll = all != nullptr;
	std::string filename = coretools::str::toString(outName, "_", header, "Histogram.txt");
	logfile().listFlush("Writing ", header, " histogram to '", filename, "' ...");
	coretools::TOutputFile out(filename, {"readGroup", header, "count"});

	TCountDistributionVector<> distributionPerReadGroup;
	coretools::TCountDistribution<> temp;
	for (size_t i = 0; i < distVec.size(); i++) {
		distVec[i].fillCombinedDistribution(temp);
		distributionPerReadGroup.add(i, temp);
	}

	if (hasAll) {
		all->fillCombinedDistribution(temp);
	}

	// Should file contain read groups with 0 counts?
	if (parameters().exists("writeZeroCounts")) {
		if (hasAll) {
			temp.write<true>(out, "allReadGroups");
		} else {
			distributionPerReadGroup.writeCombined<true>(out, "allReadGroups");
		}
		distributionPerReadGroup.write<true>(out, ReadGroupNames);
	} else {
		if (hasAll) {
			temp.write<false>(out, "allReadGroups");
		} else {
			distributionPerReadGroup.writeCombined<false>(out, "allReadGroups");
		}
		distributionPerReadGroup.write<false>(out, ReadGroupNames);
	}

	out.close();
	logfile().done();
}

size_t countsOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec) {
	size_t counts = 0;
	for (const auto &s : vec) { counts += s.counts(); }
	return counts;
}

size_t countsForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID) {
	size_t counts = 0;
	for (const auto &s : vec) { counts += s[chromRefID].counts(); }
	return counts;
}

size_t sumOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec) {
	size_t sum = 0;
	for (const auto &s : vec) { sum += s.sum(); }
	return sum;
}

size_t sumForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID) {
	size_t sum = 0;
	for (const auto &s : vec) { sum += s[chromRefID].sum(); }
	return sum;
}

double meanOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec) {
	size_t counts = countsOverAllReadGroups(vec);
	if (counts == 0) { return 0.0; }
	size_t sum = sumOverAllReadGroups(vec);
	return ((double)sum / (double)counts);
}

double meanForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID) {
	size_t counts = countsForChromosome(vec, chromRefID);
	if (counts == 0) { return 0.0; }
	size_t sum = sumForChromosome(vec, chromRefID);

	return ((double)sum / (double)counts);
}

size_t maxOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec) {
	size_t max = 0;
	for (const auto &s : vec) {
		size_t curMax = s.max();
		if (curMax > max) max = curMax;
	}
	return max;
}

size_t maxForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID) {
	size_t max = 0;
	for (const auto &s : vec) {
		size_t curMax = s[chromRefID].max();
		if (curMax > max) max = curMax;
	}
	return max;
}

size_t countsLargerZeroOverAllReadGroups(const std::vector<coretools::TCountDistributionVector<>> &vec) {
	size_t counts = 0;
	for (const auto &s : vec) { counts += s.countsLargerZero(); }
	return counts;
}

size_t countsLargerZeroForChromosome(const std::vector<coretools::TCountDistributionVector<>> &vec, size_t chromRefID) {
	size_t counts = 0;
	for (const auto &s : vec) { counts += s[chromRefID].countsLargerZero(); }
	return counts;
}

size_t count(const std::vector<std::array<size_t, 2>> &Vals, bool Paired) {
	size_t counts = 0;
	for (const auto &v : Vals) { counts += v[Paired]; }
	return counts;
}

size_t count(const std::vector<std::vector<std::array<size_t, 2>>> &Vals, bool Paired) {
	size_t counts = 0;
	for (const auto &vv : Vals) {
		for (const auto &v : vv) { counts += v[Paired]; }
	}
	return counts;
}

size_t countChr(const std::vector<std::vector<std::array<size_t, 2>>> &Vals, size_t RefID, bool Paired) {
	size_t counts = 0;
	for (const auto &vv : Vals) {
		 counts += vv[RefID][Paired]; 
	}
	return counts;
}

} // namespace impl

void TBamDiagnoser::_handleAlignment() {
	// get read group
	const auto &read    = _genome.bamFile().curRead();
	const size_t readGroup = read.readGroupID;
	if (readGroup == BAM::TReadGroups::noReadGroupId) return;

	const size_t refID = read.position.refID();

	// increments for each read that passed filters
	_passedQC.add(readGroup, refID);

	// add to counters
	constexpr size_t maxReadDist = 10000;
	const auto &curPosition = read.position;
	if (curPosition.refID() == _old[readGroup].position.refID()) {
		_readDist[readGroup].add(refID, std::min(maxReadDist, read.position - _old[readGroup].position));
		if (_identifyDuplicates && (curPosition.position() == _old[readGroup].position.position()) && (read.fragmentLength() == _old[readGroup].length)) {
			_duplicateFile.writeln(_genome.bamFile().curChromosome().name(), curPosition.position(), read.name(),
								 read.fragmentLength(), read.isReverseStrand(), _old[readGroup].name, _old[readGroup].length, _old[readGroup].isReversed);
		}
	}
	if (curPosition == _old[readGroup].position) {
		++_startCounter[readGroup];
	} else {
		_readStart[readGroup].add(refID, _startCounter[readGroup]);
		_startCounter[readGroup] = 1;
	}


	_old[readGroup].position = curPosition;
	if (_identifyDuplicates) {
		_old[readGroup].name       = read.name();
		_old[readGroup].length     = read.fragmentLength();
		_old[readGroup].isReversed = read.isReverseStrand();
	}

	if (curPosition.refID() == _oldPosition.refID()) {
		_allReadDist.add(refID, std::min(maxReadDist, read.position - _oldPosition));
	}
	if (curPosition == _oldPosition) {
		++_allStart;
	} else {
		_allReadStart.add(refID, _allStart);
		_allStart = 1;
	}
	_oldPosition = curPosition;

	_readLength[LengthType::All][readGroup].add(refID, read.cigar.lengthRead());
	if (_writeMates) {
		const auto mate1 = read.isFirstMate();
		const auto rev   = read.isReverseStrand();
		if (mate1) {
			if (rev) _readLength[LengthType::Rev1][readGroup].add(refID, read.cigar.lengthRead());
			else _readLength[LengthType::Fwd1][readGroup].add(refID, read.cigar.lengthRead());
		} else {
			if (rev) _readLength[LengthType::Rev2][readGroup].add(refID, read.cigar.lengthRead());
			else _readLength[LengthType::Fwd2][readGroup].add(refID, read.cigar.lengthRead());
		}
	}
	_usableLength[LengthType::All][readGroup].add(refID, read.cigar.lengthAligned());
	if (_writeMates) {
		const auto mate1 = read.isFirstMate();
		const auto rev   = read.isReverseStrand();
		if (mate1) {
			if (rev) _usableLength[LengthType::Rev1][readGroup].add(refID, read.cigar.lengthAligned());
			else _usableLength[LengthType::Fwd1][readGroup].add(refID, read.cigar.lengthAligned());
		} else {
			if (rev) _usableLength[LengthType::Rev2][readGroup].add(refID, read.cigar.lengthAligned());
			else _usableLength[LengthType::Fwd2][readGroup].add(refID, read.cigar.lengthAligned());
		}
	}

	_softClippedLength[readGroup].add(refID, read.cigar.lengthSoftClipped());
	_mappingQuality[readGroup].add(refID, read.mappingQuality());
	++_paired[readGroup][refID][read.isPaired()];

	// fragment length: only for proper pairs and only once
	if (read.isProperPair() && !read.isReverseStrand())
		_fragmentLength[readGroup].add(refID, read.fragmentLength());
}

	TBamDiagnoser::TBamDiagnoser() : _identifyDuplicates(parameters().exists("identifyDuplicates")), _writeMates(parameters().exists("writeMates")) {
	_genome.bamFile().readGroups().fillVectorWithNames(_readGroupNames);
}

void TBamDiagnoser::run() {
	// calculate length of genome
	if (_identifyDuplicates)
		_duplicateFile.open(_genome.outputName() + "_potentialDuplicates.txt.gz",
							{"Chr", "Pos", "Read1", "Length1", "isReversed1", "Read2", "Length2", "isReversed2"});

	// initialize counters
	const auto totLengthOfGenome = _genome.bamFile().chromosomes().referenceLength();
	const size_t numRG           = _genome.bamFile().readGroups().size();
	const size_t numChrom        = _genome.bamFile().chromosomes().size();

	// resize distributions
	_passedQC.resize(numRG);
	_passedQC.resizeDistributions(numChrom);

	_old.resize(numRG);
	_startCounter.resize(numRG, 1);
	_readDist.resize(numRG);
	_softClippedLength.resize(numRG);
	_mappingQuality.resize(numRG);
	_fragmentLength.resize(numRG);
	_readStart.resize(numRG);
	_paired.resize(numRG);

	for (auto i = LengthType::min; i < LengthType::max; ++i) {
		_readLength[i].resize(numRG);
		for (auto& rl: _readLength[i]) rl.resize(numChrom);

		_usableLength[i].resize(numRG);
		for (auto& ul: _usableLength[i]) ul.resize(numChrom);
		if (!_writeMates) break;
	}
	for (size_t i = 0; i < numRG; i++) {
		_readDist[i].resize(numChrom);
		_softClippedLength[i].resize(numChrom);
		_mappingQuality[i].resize(numChrom);
		_fragmentLength[i].resize(numChrom);
		_readStart[i].resize(numChrom);
		_paired[i].resize(numChrom);
	}
	_allReadDist.resize(numChrom);
	_allReadStart.resize(numChrom);

	if (parameters().exists("perChromosome")) {
		_chromStats = true;
		logfile().list("Will output data per chromosome into diagnostics-file. (parameter 'perChromosome')");
	} else {
		logfile().list("Will not output data per chromosome into diagnostics-file. (parameter 'perChromosome')");
	}

	// now parse through bam file
	_traverseBAMPassedQC();

	// need to add positions on the chromosome without a start
	for (const auto& chr: _genome.bamFile().chromosomes()) {
		const auto rID = chr.refID();
		for (size_t i = 0; i < numRG; ++i) {
			_readStart[i].add(rID, 0, chr.length() - _readStart[i][rID].counts());
		}
		_allReadStart.add(rID, 0, chr.length() - _allReadStart[rID].counts());
	}

	if (!parameters().exists("mergeInput")) {
		logfile().list("Will not create input file for mergeOverlappingReads. (use 'mergeInput' to do so).");
	}
	if (!parameters().exists("printReferenceLength")) {
		logfile().list(
			"Will not print reference lengths of chromosomes to file. (use 'printReferenceLength' to do so).");
	}
	logfile().list("Approximate sequencing depth was estimated at ",
				   (double)impl::sumOverAllReadGroups(_usableLength[LengthType::All]) / totLengthOfGenome, ".");

	// writing output files
	logfile().startIndent("Writing output files:");

	// writing read group summary
	std::string filename = _genome.outputName() + "_diagnostics.txt";
	logfile().listFlush("Writing general diagnostics to '" + filename + "' ...");
	coretools::TOutputFile out(filename);

	std::vector<std::string_view> header{"readGroup"};
	if (_chromStats) header.push_back("chromosome");
	header.insert(header.end(),
				  {"totalReads", "passedQC", "duplicates", "avgReadLength", "seqCycles", "avgReadDist", "avgReadStart",
				   "singleEnd", "pairedEnd", "properPairs", "avgFragmentLength", "softClipped", "avgSoftClippedLength",
				   "avgUsableAlignedLength", "approximateDepth", "avgMappingQuality", "seqType"});
	out.writeHeader(header);

	// write for all read groups and all chromosomes
	out.write("allReadGroups");
	const auto singlesAll = impl::count(_paired, false);
	const auto pairsAll   = impl::count(_paired, true);
	if (_chromStats) { out.write("allChromosomes"); }
	out.write(_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome().counts(), _passedQC.counts(),
			  _genome.bamFile().filter(BAM::FilterType::Duplicate).getCombinedCounts(),
			  impl::meanOverAllReadGroups(_readLength[LengthType::All]), impl::maxOverAllReadGroups(_readLength[LengthType::All]), _allReadDist.mean(),
			  _allReadStart.mean(), singlesAll, pairsAll,
			  impl::countsOverAllReadGroups(_fragmentLength), impl::meanOverAllReadGroups(_fragmentLength),
			  impl::countsLargerZeroOverAllReadGroups(_softClippedLength),
			  impl::meanOverAllReadGroups(_softClippedLength), impl::meanOverAllReadGroups(_usableLength[LengthType::All]),
			  (double)impl::sumOverAllReadGroups(_usableLength[LengthType::All]) / totLengthOfGenome,
			  impl::meanOverAllReadGroups(_mappingQuality));
	if (singlesAll > 0 && pairsAll > 0) {
		out.writeln("mixed");
	} else if (singlesAll > 0) {
		out.writeln("single");
	} else if (pairsAll > 0) {
		out.writeln("paired");
	} else {
		out.writeln("empty");
	}

	// write for all read groups per chromosome
	if (_chromStats) {
		for (const auto &chr : _genome.bamFile().chromosomes()) {
			size_t refID       = chr.refID();
			const auto singles = impl::countChr(_paired, refID, false);
			const auto pairs   = impl::countChr(_paired, refID, true);
			out.write("allReadGroups", chr.name(),
					  (_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome()).horizontalCounts(refID),
					  _passedQC.horizontalCounts(refID),
					  _genome.bamFile().filter(BAM::FilterType::Duplicate).getCountsPerChromosome(refID),
					  impl::meanForChromosome(_readLength[LengthType::All], refID), impl::maxForChromosome(_readLength[LengthType::All], refID),
					  _allReadDist[refID].mean(), _allReadStart[refID].mean(), singles, pairs,
					  impl::countsForChromosome(_fragmentLength, refID),
					  impl::meanForChromosome(_fragmentLength, refID),
					  impl::countsLargerZeroForChromosome(_softClippedLength, refID),
					  impl::meanForChromosome(_softClippedLength, refID), impl::meanForChromosome(_usableLength[LengthType::All], refID),
					  (double)impl::sumForChromosome(_usableLength[LengthType::All], refID) / (double)chr.length(),
					  impl::meanForChromosome(_mappingQuality, refID));
			if (singles > 0 && pairs > 0) {
				out.writeln("mixed");
			} else if (singles > 0) {
				out.writeln("single");
			} else if (pairs > 0) {
				out.writeln("paired");
			} else {
				out.writeln("empty");
			}
		}
	}

	// write per read group for all chromosomes
	for (size_t rg = 0; rg < numRG; ++rg) {
		const auto singles = impl::count(_paired[rg], false);
		const auto pairs   = impl::count(_paired[rg], true);

		out.write(_genome.bamFile().readGroups().getName(rg));
		if (_chromStats) out.write("allChromosomes");
		out.write((_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome())[rg].counts(), _passedQC[rg].counts(),
		          _genome.bamFile().filter(BAM::FilterType::Duplicate).getCounts(rg), _readLength[LengthType::All][rg].mean(),
				  _readLength[LengthType::All][rg].max(), _readDist[rg].mean(), _readStart[rg].mean(), singles, pairs,
		          _fragmentLength[rg].counts(), _fragmentLength[rg].mean(), _softClippedLength[rg].countsLargerZero(),
		          _softClippedLength[rg].mean(), _usableLength[LengthType::All][rg].mean(),
		          (double)_usableLength[LengthType::All][rg].sum() / (double)totLengthOfGenome, _mappingQuality[rg].mean());

		if (singles > 0 && pairs > 0) {
			out.writeln("mixed");
		} else if (singles > 0) {
			out.writeln("single");
		} else if (pairs > 0) {
			out.writeln("paired");
		} else {
			out.writeln("empty");
		}
		if (_chromStats) {
			// write per read group per chromosome
			for (const auto& chr: _genome.bamFile().chromosomes()) {
				size_t refID       = chr.refID();
				const auto singles = _paired[rg][refID][false];
				const auto pairs   = _paired[rg][refID][true];
				out.write(
				    _genome.bamFile().readGroups().getName(rg), chr.name(),
				    (_genome.bamFile().numAlignmentReadPerReadGroupPerChromosome())[rg][refID], _passedQC[rg][refID],
				    _genome.bamFile().filter(BAM::FilterType::Duplicate).getCountsAtReadGroupAndChromosome(rg, refID),
				    _readLength[LengthType::All][rg][refID].mean(), _readLength[LengthType::All][rg][refID].max(), _readDist[rg][refID].mean(),
				    _readStart[rg][refID].mean(), singles, pairs,
				    _fragmentLength[rg][refID].counts(), _fragmentLength[rg][refID].mean(),
				    _softClippedLength[rg][refID].countsLargerZero(), _softClippedLength[rg][refID].mean(),
				    _usableLength[LengthType::All][rg][refID].mean(), (double)_usableLength[LengthType::All][rg][refID].sum() / (double)chr.length(),
				    _mappingQuality[rg][refID].mean());

				if (singles > 0 && pairs > 0) {
					out.writeln("mixed");
				} else if (singles > 0) {
					out.writeln("single");
				} else if (pairs > 0) {
					out.writeln("paired");
				} else {
					out.writeln("empty");
				}
			}
		}
	}
	out.close();
	logfile().done();

	if (parameters().exists("printReferenceLength")) {
		// write file with length of all contigs
		std::string referenceLengthName = _genome.outputName() + "_referenceLengths.txt";
		logfile().listFlush("Outputting reference lengths of all contigs to '" + referenceLengthName + "' ...");
		coretools::TOutputFile refLen(referenceLengthName, {"chromosome", "length"});
		auto it = _genome.bamFile().chromosomes().cbegin();
		while (it != _genome.bamFile().chromosomes().cend()) {
			refLen << it->name() << it->length() << coretools::endl;
			++it;
		}
		refLen.close();
		logfile().done();
	}

	// writing distributions
	impl::writeHistogram(_readLength[LengthType::All], _genome.outputName(), "readLength", _readGroupNames);
	if (_writeMates) {
		impl::writeHistogram(_readLength[LengthType::Fwd1], _genome.outputName(), "readLengthFwd1", _readGroupNames);
		impl::writeHistogram(_readLength[LengthType::Fwd2], _genome.outputName(), "readLengthFwd2", _readGroupNames);
		impl::writeHistogram(_readLength[LengthType::Rev1], _genome.outputName(), "readLengthRev1", _readGroupNames);
		impl::writeHistogram(_readLength[LengthType::Rev2], _genome.outputName(), "readLengthRev2", _readGroupNames);
	}
	impl::writeHistogram(_readDist, _genome.outputName(), "readDist", _readGroupNames, &_allReadDist);
	impl::writeHistogram(_readStart, _genome.outputName(), "readStart", _readGroupNames, &_allReadStart);
	impl::writeHistogram(_usableLength[LengthType::All], _genome.outputName(), "alignedLength", _readGroupNames);
	if (_writeMates) {
		impl::writeHistogram(_usableLength[LengthType::Fwd1], _genome.outputName(), "alignedLengthFwd1", _readGroupNames);
		impl::writeHistogram(_usableLength[LengthType::Fwd2], _genome.outputName(), "alignedLengthFwd2", _readGroupNames);
		impl::writeHistogram(_usableLength[LengthType::Rev1], _genome.outputName(), "alignedLengthRev1", _readGroupNames);
		impl::writeHistogram(_usableLength[LengthType::Rev2], _genome.outputName(), "alignedLengthRev2", _readGroupNames);
	}
	impl::writeHistogram(_softClippedLength, _genome.outputName(), "softClippedLength", _readGroupNames);
	impl::writeHistogram(_fragmentLength, _genome.outputName(), "fragmentLength", _readGroupNames);
	impl::writeHistogram(_mappingQuality, _genome.outputName(), "mappingQuality", _readGroupNames);

	logfile().endIndent(); // end writing output files
}

} // namespace GenomeTasks
