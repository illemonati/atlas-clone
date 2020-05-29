/*
 * TAlignmentMerger.cpp
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#include "TAlignmentMerger.h"

//-----------------------------------------
// TMateFinder
//-----------------------------------------
TMateFilter::TMateFilter(BAM::TBamFile& BamFile, TParameters & params, TLog* Logfile):_bamFile(BamFile){
	logfile = Logfile;

	//max distance between mates
	_maxDistanceBetweenMates = params.getParameterIntWithDefault("acceptedDistance", 2000);
	logfile->list("Mates that are farther than " + toString(_maxDistanceBetweenMates) + " apart will be considered orphans. (parameter 'acceptedDistance')");

	//keep orphans
	if(params.parameterExists("keepOrphans")){
		_keepOrphans = true;
		logfile->list("Will keep keep orphaned reads. (parameter 'keepOrphans')");
	} else {
		_keepOrphans = false;
		logfile->list("Will filter out orphaned reads (use 'keepOrphans' to keep them).");
	}
};

TAlignmentInStorage TMateFilter::_findMate(const std::string & name){
	for(auto it=_alignmentStorage.begin(); it!=_alignmentStorage.end(); ++it){
		if(it->alignment->_name() == name){
			return it;
		}
	}

	return _alignmentStorage.end();
};


void TMateFilter::_writeAlignment(TAlignmentInStorage & it){
	//save the alignment to the bam file
	_bamFile.writeAlignment(*(it->alignment), genoMap, qualMap);
	//delete it->alignment;
	it = _alignmentStorage.erase(it);
};

void TMateFilter::_writeOrFilterAsOrphan(TAlignmentInStorage & it){
	if(it->ready){
		_writeAlignment(it);
	} else if(_keepOrphans){
		//set as improper pair
		it->setAsNonProperPair();
		//write to BAM file
		_writeAlignment(it);
	} else {
		//write reason to bam log
		_bamFile.filterOut(it->alignment->name(), it->alignment->isReverseStrand());
		//delete it->alignment;
		it = _alignmentStorage.erase(it);
	}
};

void TMateFilter::_writeAll(){
	//write everything and mark reads with missing mates as improper.
	//reads still in storage are no-proper pairs: write or add to black list
	TAlignmentInStorage it = _alignmentStorage.begin();
	while(it != _alignmentStorage.end()){
		_writeOrFilterAsOrphan(it);
	}
	//clear blacklist: future reads will anyways be orphans
	_blacklist.clear();
};

void TMateFilter::_writeUpTo(const int position){
	//writes all that are ready or too far away
	TAlignmentInStorage it = _alignmentStorage.begin();
	while(it != _alignmentStorage.end() && (it->ready || position - it->alignment->position() > _maxDistanceBetweenMates)){
		_writeOrFilterAsOrphan(it);
	}
};

void TMateFilter::_handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate){
	if(!alignment->isProperPair()){
		//not a proper pair: mark mate as as improper
		mate->setAsNonProperPair();
	}
	//mark both as ready for writing
	mate->ready = true;
	_alignmentStorage.emplace_back(alignment, true);
};

void TMateFilter::traverseBAM(const std::string outputName){
	//open writer
	_bamFile.openOutput(outputName + "_filtered.bam");
	_bamFile.setExternalFilterReason("Orphan");

	//now parse BAM file
	_bamFile.startProgressReporting(1000000, logfile);
	while(_bamFile.readNextAlignment()){
		//if on new chromosome, empty storage
		if(_bamFile.chrChanged()){
			//write all ready currently in storage
			_writeAll();
		}

		//check if first alignment in storage is too far away from current alignment
		//if yes, first alignment in storage is considered an orphan
		_writeUpTo(_bamFile.curPosition());

		//check if read passed filters
		if(_bamFile.curPassedQC()){
			//if single end and storage is empty: write directly
			if(!_bamFile.curIsPaired() && _alignmentStorage.empty()){
				_bamFile.writeCurAlignment();
			} else {
				//parse alignment
				BAM::TAlignment* alignment;
				_bamFile.fill(*alignment);

				//if read is paired, check for mate
				if(alignment->isPaired()){
					//if mate is in blacklist: add as improper pair for writing
					if(_blacklist.isInBlacklist(alignment->name())){
						alignment->setIsProperPair(false);
						_alignmentStorage.emplace_back(alignment, true);
						_blacklist.remove(alignment->name());
					} else {
						//check if mate is in storage. If so, add as ready to write too
						auto mate = _findMate(alignment->name());
						if(mate == _alignmentStorage.end()){
							//no mate found
							if(alignment->isProperPair()){
								//proper pair: add to storage and wait for mate
								_alignmentStorage.emplace_back(alignment, false);
							} else {
								//improper pair: add to blacklist and ready to write
								_blacklist.add(alignment->name());
								_alignmentStorage.emplace_back(alignment, true);
							}
						} else {
							//mate found
							if(!alignment->isProperPair()){
								//not a proper pair: mark mate as as improper
								mate->setAsNonProperPair();
							}
							//mark both as ready for writing
							mate->ready = true;
							_alignmentStorage.emplace_back(alignment, true);
						}
					}
				} else {
					//read is single end: add for writing
					_alignmentStorage.emplace_back(alignment, true);
				}
			}
		} else {
			//Did not pass QC: filter out
			//need to store in blacklist if is was paired
			if(_bamFile.curIsProperPair()){
				_blacklist.add(_bamFile.curName());
			}
		}

		//report
		_bamFile.printProgress();
	}

	//write reads still in storage
	_writeAll();

	//done parsing bam file: report
	_bamFile.printSummary(logfile);
};

//-----------------------------------------
// TAlignmentMergerType
//-----------------------------------------
TAlignmentMergerType_none::TAlignmentMergerType_none(){
	_numOverlap = 0;
};

void TAlignmentMergerType_none::_prepare(){
	_numOverlap = 0;
};

void TAlignmentMergerType_none::_mergeBases(TBase & bestBase, TBase & worstBase){
	if(_adaptQuality){
		GenotypeLikelihoods::TBaseData likelihood(bestBase.base, qualMap.phredIntToError(bestBase.recalibratedQualityAsPhredInt));
		likelihood *= GenotypeLikelihoods::TBaseData(worstBase.base, qualMap.phredIntToError(worstBase.recalibratedQualityAsPhredInt));
		likelihood.normalize();
		bestBase.recalibratedQualityAsPhredInt = qualMap.errorToPhredInt(1.0 - likelihood[bestBase.base]);
	}

	//set other to missing
	worstBase.recalibratedQualityAsPhredInt = 0.0;
	worstBase.base = N;
};

void TAlignmentMergerType_none::merge(BAM::TAlignment* alignment, BAM::TAlignment* mate){
	//NOTE: mate is earlier!
	//deletions and insertions are kept as is. these positions are not compared

	//prepare (e.g. pcik random number
	_prepare();

	//go through alignments
	int fwdP = 0; int revP = 0;
	while(fwdP <= alignment->lastAlingedInternalPos() && revP <= mate->lastAlingedInternalPos()){
		//make sure we comapre at the same position in respect to ref
		if(!alignment->isAlignedAtInternalPos(fwdP)){
			++fwdP;
		} else if(!mate->isAlignedAtInternalPos(revP)){
			++revP;
		} else if(alignment->positionInRef(fwdP) < mate->positionInRef(revP)){
			++fwdP;
		} else if(alignment->positionInRef(fwdP) > mate->positionInRef(revP)){
			++revP;
		} else {
			//at same position: merge
			if(_merge()){
				++_numOverlap;
			}
			++fwdP; ++revP;
		}


		if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
			//bases overlap same position in ref -> choose at random which to keep
			if(keepFwd){
				_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
			} else {
				_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
			}
			//increment both counters
			++fwdP;
			++revP;
		} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
			++fwdP;
		} else {
			++revP;
		}
	}
};

TAlignmentMergerType_randomRead::TAlignmentMergerType_randomRead(TRandomGenerator* RandomGenerator){
	_randomGenerator = RandomGenerator;
	_keepMate = false;
};

void TAlignmentMergerType_randomRead::_prepare(){
	_numOverlap = 0;
	_keepMate = _randomGenerator->pickOneOfTwo();
};


bool TAlignmentMergerType_randomRead::_merge(BAM::TAlignment* alignment, BAM::TAlignment* mate){

};

//-----------------------------------------
// TAlignmentMerger
//-----------------------------------------
TAlignmentMerger::TAlignmentMerger(BAM::TBamFile& BamFile, TParameters & params, TLog* Logfile, TRandomGenerator* RandomGenerator):TMateFilter(BamFile, params, Logfile){
	_randomGenerator = RandomGenerator;


	//check if keepAllReads is turned on
	//TODO: what is the basic set of filters needed?
	if(!_bamFile.improperPairsFilterEnabled()){
		logfile->warning("Improper pairs are kept but will not be merged!");
	}
	/*
	if(alignmentParser.getKeepAll()){
		logfile->warning("Undefined behavior when merging reads that do not pass default filters. Consider removing 'keepAllReads'");
	}
	*/

	//set merging method
	//TODO: update wiki to relfect change in names
	std::string method = params.getParameterStringWithDefault("mergingMethod", "keepRandomRead");
	if(method == "none"){
		_type = none;
		logfile->list("Merging method: no merging.");
	} else if (method == "randomRead"){
		_type = randomRead;
		logfile->list("Merging method: will keep random read for all overlapping positions");
	} else if(method == "randomBase"){
		_type = randomBase;
		logfile->list("Merging method: will keep random base at each overlapping positios.");
	} else if(method == "highestQuality"){
		_type = highestQuality;
		logfile->list("Merging method: will keep base with highest quality at overlapping positions.");
	} else {
		throw "Unknown merging method " + method + "! Use 'none', 'randomRead', 'randomBase' and 'highestQuality'.";
	}

	//decide if we update quality score
	if(params.parameterExists("updateQuality")){
		_adaptQuality = true;
		logfile->list("Will update quality scores of prefered bases to reflect information from overlapping bases.");
	} else {
		_adaptQuality = false;
		logfile->list("Will keep original quality scores of the prefered bases (use updateQuality to update quality scores).");
	}

};

void TAlignmentMerger::_handleMates(BAM::TAlignment* alignment, TAlignmentInStorage & mate){
	if(!alignment->isProperPair()){
		//not a proper pair: mark mate as as improper too
		mate->setAsNonProperPair();
	} else {
		//attempt merging
		if(_type != none && mate->alignment->lastAlignedPositionWithRespectToRef() >= alignment->position()){
			//reads overlap: merge
			switch(_type){
			case randomRead:

				break;
			case randomBase:

				break;
			case highestQuality:

				break;
			}

			//set alignments have changes
			alignment->setHasChanged();
			mate->alignment->setHasChanged();
		}
	}

	//mark both as ready for writing
	mate->ready = true;
	_alignmentStorage.emplace_back(alignment, true);
};




void TAlignmentMerger::_mergeRandomRead(BAM::TAlignment* alignment, BAM::TAlignment* mate){
	//NOTE: mate is earlier!
	//deletions and insertions are kept as is. these positions are not compared

	//which read to keep
	bool keepFwd = _randomGenerator->pickOneOfTwo();

	//reads overlap -> check if there are bases overlapping same position in ref
	//alignedPos is with respect to read
	int fwdP = 0; int revP = 0;
	while(fwdP <= alignment->lastAlingedInternalPos() && revP <= mate->lastAlingedInternalPos()){



		if(!alignment->isAlignedAtInternalPos(fwdP)){
			++fwdP;
		} else if(!mate->isAlignedAtInternalPos(revP)){
			++revP;
		} else if(alignment->positionInRef(fwdP) < mate->positionInRef(revP)){
			++fwdP;
		} else if(alignment->positionInRef(fwdP) > mate->positionInRef(revP)){
			++revP;
		} else {
			//at same position: merge

		}


		if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
			//bases overlap same position in ref -> choose at random which to keep
			if(keepFwd){
				_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
			} else {
				_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
			}
			//increment both counters
			++fwdP;
			++revP;
		} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
			++fwdP;
		} else {
			++revP;
		}
	}
}
}

void TAlignmentParser::mergeAlignedBasesBamReadsRandom(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality, TRandomGenerator* randomGenerator){
	//deletions and insertions are kept as is. these positions are not compared
	if(fwdAlignment->lastAlignedPositionWithRespectToRef >= revAlignment->position){
		fwdAlignment->setAlignmentHasChanged();
		revAlignment->setAlignmentHasChanged();

		//reads overlap -> check if there are bases overlapping same position in ref
		//alignedPos is with respect to read
		int fwdP = 0;
		int revP = 0;
		while(fwdP <= fwdAlignment->lastAlignedPos && revP <= revAlignment->lastAlignedPos){
			if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
				//bases overlap same position in ref -> choose at random which to keep
				if(randomGenerator->getRand() < 0.5){
					_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
};

void TAlignmentParser::mergeAlignedBasesBamReads(TAlignment* fwdAlignment, TAlignment* revAlignment, bool adaptQuality){
	//deletions and insertions are kept as is. these positions are not compared
	if(fwdAlignment->lastAlignedPositionWithRespectToRef >= revAlignment->position){
		fwdAlignment->setAlignmentHasChanged();
		revAlignment->setAlignmentHasChanged();

		//reads overlap -> check if there are bases overlapping same position in ref
		//alignedPos is with respect to read
		int fwdP = 0;
		int revP = 0;
		while(fwdP <= fwdAlignment->lastAlignedPos && revP <= revAlignment->lastAlignedPos){
			if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] == revAlignment->position + revAlignment->alignedPosition[revP]){
				//bases overlap same position in ref -> keep the one with higher quality
				if(fwdAlignment->bases[fwdP].recalibratedQualityAsPhredInt > revAlignment->bases[revP].recalibratedQualityAsPhredInt){
					_adaptQualityWhenMerging(fwdAlignment->bases[fwdP], revAlignment->bases[revP], adaptQuality);
				} else {
					_adaptQualityWhenMerging(revAlignment->bases[revP], fwdAlignment->bases[fwdP], adaptQuality);
				}
				//increment both counters
				++fwdP;
				++revP;
			} else if(fwdAlignment->position + fwdAlignment->alignedPosition[fwdP] < revAlignment->position + revAlignment->alignedPosition[revP]){
				++fwdP;
			} else {
				++revP;
			}
		}
	}
};


void TMateFilter::mergeBAM(const std::string outputName){
	_traverseBAM(outputName);
};


TAlignmentMerger::TAlignmentMerger(BamTools::BamWriter* Writer, TAlignmentParser* Parser, int MaxDistanceBetweenMates){
	writer = Writer;
	parser = Parser;
	_maxDistanceBetweenMates = MaxDistanceBetweenMates;
	_keepHigherQuality = false;
	_adaptQuality = false;
	_filterOrphans = true;
	_keepRandomBase = false;
	_keepRandomRead = false;
	_allowForLarger = false;
};

void TAlignmentMerger::_writeAlignment(std::vector< TAlignmentMergerEntry >::iterator & it){
	//save the alignment to the bam file
	bamFile.writeAlignment(*(it->alignment), genoMap, qualMap);
	//delete it->alignment;
	it = alignmentStorage.erase(it);
};

void TAlignmentMerger::_addToBlacklist(std::vector< TAlignmentMergerEntry >::iterator & it, std::string error){
	parser->addToBlacklist(*(it->alignment), error);
	//delete it->alignment;
	it = alignmentStorage.erase(it);
};

void TAlignmentMerger::_writeAllThatAreReady(){
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();
	while(it != alignmentStorage.end() && it->ready){
		_writeAlignment(it);
	}
};

std::vector< TAlignmentMergerEntry >::iterator TAlignmentMerger::_findMate(TAlignment & alignment){
	std::vector< TAlignmentMergerEntry >::iterator it;
	for(it=alignmentStorage.begin(); it!=alignmentStorage.end(); ++it){
		//found its mate!
		if(it->alignment->_name() == alignment._name()){
			return it;
		}
	}

	return alignmentStorage.end();
};

void TAlignmentMerger::addToBeMerged(TAlignment & alignment, TRandomGenerator* randomGenerator){
	std::vector< TAlignmentMergerEntry >::iterator it = _findMate(alignment);
	if(it == alignmentStorage.end()){
		//no mate found: add to storage
		alignmentStorage.emplace_back(alignment, false);
	} else {
		//mate found, merge!
		if(_keepRandomRead)
			parser->mergeAlignedBasesOneRead(it->alignment, &alignment, _adaptQuality, randomGenerator);
		if(_keepRandomBase)
			parser->mergeAlignedBasesBamReadsRandom(it->alignment, &alignment, _adaptQuality, randomGenerator);
		else
			parser->mergeAlignedBasesBamReads(it->alignment, &alignment, _adaptQuality);
		it->ready = true;
		alignmentStorage.emplace_back(alignment, true);
	}
};

void TAlignmentMerger::addToBeSplit(TAlignment & alignment, std::map<int, TReadGroupMaxLength>::iterator singleEndRGIT){
	//check length
	if(alignment.getBamAlignmentLength() == singleEndRGIT->second.maxLen){
		alignment.updateOptionalSamField("RG", singleEndRGIT->second.truncatedOrMergedReadGroup);
	} else if(alignment.getBamAlignmentLength() > singleEndRGIT->second.maxLen){
		if(_allowForLarger)
			alignment.updateOptionalSamField("RG", singleEndRGIT->second.truncatedOrMergedReadGroup);
		else {
			throw("Length of read " + alignment.name() + " is > max length provided for its read group (" + toString(singleEndRGIT->second.maxLen) + ")! Use parameter 'allowForLarger' to ignore and put read in truncated read group.");
		}
	}

	//add to storage
	alignmentStorage.emplace_back(alignment, true);
}

void TAlignmentMerger::checkForMateAndWriteUnmerged(TAlignment & alignment){
	std::vector< TAlignmentMergerEntry >::iterator it = _findMate(alignment);
	if(it == alignmentStorage.end()){
		//no mate found: add to storage
		alignmentStorage.emplace_back(alignment, false);
	} else {
		//mate found, ready to write!
		it->ready = true;
		alignmentStorage.emplace_back(alignment, true);
	}
};

void TAlignmentMerger::addReadyToBeWritten(TAlignment & alignment){
	if(alignmentStorage.empty()){
		bamFile.writeAlignment(alignment, genoMap, qualMap);
	} else {
		alignmentStorage.emplace_back(alignment, true);
	}
};

void TAlignmentMerger::addAsImproperPair(TAlignment & alignment){
	if(_filterOrphans){
		//no need to keep mate in list anymore
		parser->removeFromBlacklist(alignment, "not a proper pair (orphan)");
	} else {
		//set to improper read
		alignment.setIsProperPair(false);
		addReadyToBeWritten(alignment);
	}
};

void TAlignmentMerger::writeUpTo(const int position){
	//writes all that are ready or too far away
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();
	while(it != alignmentStorage.end() && (it->ready || position - it->alignment->getPosition() > _maxDistanceBetweenMates)){
//		std::cout << "it->alignment->name " << it->alignment->alignmentName << " position - it->alignment->position " << position << "-" <<  it->alignment->position << std::endl;
		if(it->ready){
			_writeAlignment(it);
		} else {
			if(_filterOrphans){
				_addToBlacklist(it, "orphaned read: mate is farther away than " + toString(_maxDistanceBetweenMates) + " bp");
			} else {
				it->setAsNonProperPair();
				_writeAlignment(it);
			}
		}
	}
};

void TAlignmentMerger::clear(){
	//write everything and mark reads with missing mates as improper.
	std::vector< TAlignmentMergerEntry >::iterator it = alignmentStorage.begin();

	//reads still in storage are no-proper pairs: write or add to black list
	while(it != alignmentStorage.end()){
		if(it->ready){
			_writeAlignment(it);
		} else {
			if(_filterOrphans){
				_addToBlacklist(it, "mate on different chromosome");
			} else {
				//set reads in storage to improper pairs but ready for writing
				it->setAsNonProperPair();
				_writeAlignment(it);
			}
		}
	}
};
