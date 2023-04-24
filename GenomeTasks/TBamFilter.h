/*
 * TAlignmentMerger.h
 *
 *  Created on: Mar 28, 2019
 *      Author: wegmannd
 */

#ifndef TBAMFILTER_H_
#define TBAMFILTER_H_

#include <functional>
#include <memory>
#include <set>
#include <stdint.h>
#include <string>

#include "TBamFile.h"
#include "TGenome.h"
#include "coretools/Main/TTask.h"
#include "coretools/Strings/stringFunctions.h"

namespace BAM { class TAlignment; }
namespace BAM { class TReadGroups; }
namespace BAM { class TSequencedBase; }
namespace genometools { class TGenomePosition; }

namespace GenomeTasks{

namespace BamFilter{

using coretools::instances::parameters;
using coretools::instances::logfile;
using namespace coretools::str;

//-----------------------------------------
// TAlignmentMergerEntry
//-----------------------------------------
class TAlignmentMergerEntry{
private:
	mutable BAM::TAlignment* _alignment;
	mutable bool _ready;

public:

	TAlignmentMergerEntry(BAM::TAlignment* Alignment, bool readyForWriting);
	TAlignmentMergerEntry(TAlignmentMergerEntry && other);
	~TAlignmentMergerEntry();
	TAlignmentMergerEntry& operator=(TAlignmentMergerEntry && other);

	//getters
	bool ready() const { return _ready; }
	const BAM::TAlignment& alignment() const { return *_alignment; }
	BAM::TAlignment* stealAlignment() const { BAM::TAlignment* tmp = _alignment; _alignment = nullptr; return tmp; }
	const std::string& name() const;

	//setters
	void makeReady() const { _ready = true; }
	void setAsNonProperPair() const;
	bool operator<(const TAlignmentMergerEntry & other) const;
};

//-----------------------------------------
// TAlignmentStorage
//-----------------------------------------
typedef std::vector< TAlignmentMergerEntry > TAlignmentStorage;
typedef std::vector< TAlignmentMergerEntry >::iterator TAlignmentStorageIterator;

void addToContainer(TAlignmentStorage & Storage, BAM::TAlignment* Alignment, bool readyForWriting);

typedef std::multiset< TAlignmentMergerEntry, std::less<TAlignmentMergerEntry> > TAlignmentStorageSorted;
typedef std::multiset< TAlignmentMergerEntry, std::less<TAlignmentMergerEntry> >::iterator TAlignmentStorageSortedIterator;

void addToContainer(TAlignmentStorageSorted & Storage, BAM::TAlignment* Alignment, bool readyForWriting);

template<typename StorageType>
auto findInStorage(StorageType & Storage, const std::string & name){
	for(auto it=Storage.begin(); it!=Storage.end(); ++it){
		if(it->name() == name){
			return it;
		}
	}
	return Storage.end();
}

//-----------------------------------------
// TGenomeParsedWithAlignmentStorage
//-----------------------------------------
template<typename StorageType, typename StorageIteratorType>
class TGenomeParsedWithAlignmentStorage:public TGenome_parsed{
protected:
	BAM::TAlignmentList _blacklist; //used to keep track of filtered out mates
	StorageType _alignmentStorage;

	BAM::TOutputBamFile _outBam;

	size_t _maxDistanceBetweenMates;
	bool _recalibrate;
	bool _incorporatePMD;
	bool _keepOrphans;
	bool _removeSoftClippedBases;
	size_t _maxNumberOfSoftClippedBases;

	void _writeAlignment(StorageIteratorType & it){
		//save the alignment to the bam file
		_outBam.writeAlignment(it->alignment());
		//delete it->alignment;
		it = _alignmentStorage.erase(it);
	}

	void _writeOrFilterAsOrphan(StorageIteratorType & it){
		if(it->ready()){
			_writeAlignment(it);
		} else if(_keepOrphans){
			//set as improper pair
			it->setAsNonProperPair();
			//write to BAM file
			_writeAlignment(it);
		} else {
			//write reason to bam log
			_bamFile.filterOut(it->alignment().name(), it->alignment().isSecondMate(), it->alignment().readGroupId(), it->alignment().refID());
			it = _alignmentStorage.erase(it);
		}
	}

	void _writeAll(){
		//write everything and mark reads with missing mates as improper.
		//reads still in storage are no-proper pairs: write or add to black list
		auto it = _alignmentStorage.begin();
		while(it != _alignmentStorage.end()){
			_writeOrFilterAsOrphan(it);
		}
		//clear blacklist: future reads will anyways be orphans
		_blacklist.clear();
	}

	void _writeUpTo(const genometools::TGenomePosition & position){
		//writes all that are ready or too far away
		auto it = _alignmentStorage.begin();
		while(it != _alignmentStorage.end() &&
			  it->alignment() < position &&
			  (it->ready() || (position - it->alignment()) > _maxDistanceBetweenMates)
	    ){
			_writeOrFilterAsOrphan(it);
		}
	}

	BAM::TAlignment* _parseIntoNewAlignment(){
		BAM::TAlignment* alignment = new BAM::TAlignment;
		_bamFile.fill(*alignment);
		if(_recalibrate){
			if(_incorporatePMD){
				alignment->recalibrateWithPMD(_genotypeLikelihoodCalculator);
			} else {
				alignment->parse(_genotypeLikelihoodCalculator.sequencingErrorModels());
			}
		}
		return alignment;
	}

	//overriding _handleAligment from TGenome to do nothing
	void _handleAlignment() override {}

	//pure virtual functions
	virtual void _openBamFileForWriting() = 0;
	virtual void _handleMates(BAM::TAlignment & alignment, StorageIteratorType mate) = 0;
	virtual void _handleSingle(BAM::TAlignment & alignment) = 0;
	virtual bool _alignmentCanBeWrittenUnchanged() = 0;


public:
	TGenomeParsedWithAlignmentStorage(){
		//max distance between mates
		_maxDistanceBetweenMates = parameters().getParameterWithDefault<int>("acceptedDistance", 2000);
		logfile().list("Mates that are farther than " + toString(_maxDistanceBetweenMates) + " apart will be considered orphans. (parameter 'acceptedDistance')");

		//keep orphans
		if(parameters().parameterExists("keepOrphans")){
			_keepOrphans = true;
			logfile().list("Will keep orphaned reads. (parameter 'keepOrphans')");
		} else {
			_keepOrphans = false;
			logfile().list("Will filter out orphaned reads. (use 'keepOrphans' to keep them)");
		}

		//recalibrate BAM?
		if(_genotypeLikelihoodCalculator.recalibrationChangesQualities() || parameters().parameterExists("incorporatePMD")){
			_recalibrate = true;
			logfile().list("Will write recalibrated quality scores.");
			if(parameters().parameterExists("incorporatePMD")){
				logfile().list("Probability of PMD will be reflected in new quality scores. (parameter 'incorporatePMD')");
				_incorporatePMD = true;
				if(!_genotypeLikelihoodCalculator.hasPMD()){
					UERROR("No PMD probabilities provided! Provide PMD probabilities or remove parameter 'incorporatePMD'.");
				}
			} else {
				_incorporatePMD = false;
				logfile().list("PMD will not be reflected in the quality scores. (recommended option. Use 'incorporatePMD' to overrule)");
			}
		} else {
			logfile().list("Will write original quality scores. (provide recalibration parameters to update quality scores)");
			_recalibrate = false;
			_incorporatePMD = false;
		}

		if (parameters().parameterExists("removeSoftClippedBases")){
			_removeSoftClippedBases = true;
			//if parameter is set and a number is given -> use this as max number of softclipped bases, else remove all
			if(!parameters().getParameter<std::string>("removeSoftClippedBases").empty()){
				_maxNumberOfSoftClippedBases = parameters().getParameter<size_t>("removeSoftClippedBases");
				logfile().list("Will leave up to " + toString(_maxNumberOfSoftClippedBases) + " softclipped bases per end. (parameter 'removeSoftClippedBases')");
			} else {
				_maxNumberOfSoftClippedBases = 0;
				logfile().list("Will remove all softclipped bases. (parameter 'removeSoftClippedBases')");
			}
		} else {
			logfile().list("Will not remove softclipped bases. (Use parameter 'removeSoftClippedBases' to do so)");
			_removeSoftClippedBases = false;
		}

	}

	virtual void traverseBAM(){
		//open writer
		_openBamFileForWriting();
		_bamFile.setExternalFilterReason("Orphan");

		//now parse BAM file
		_bamFile.startProgressReporting(1000000);
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
				//if single end, unchanged and storage is empty: write directly
				if(_alignmentCanBeWrittenUnchanged()){
					_bamFile.writeCurAlignment(_outBam);
				} else {
					//parse alignment
					BAM::TAlignment* alignment = _parseIntoNewAlignment();

					if (_removeSoftClippedBases) {
						// parse and then remove softclipped reads
						alignment->parse();
						alignment->removeSoftClippedBases(_maxNumberOfSoftClippedBases);
					} 

					//if read is paired, check for mate
					if(alignment->isPaired()){
						//if mate is in blacklist: add as improper pair for writing
						if(_blacklist.isInBlacklist(alignment->name())){
							//TODO: should we mark them as improper or not?? Are all in blacklist already improper pair?
							//alignment->setIsProperPair(false);
							addToContainer(_alignmentStorage, alignment, false);
							_blacklist.remove(alignment->name());
						} else {
							//check if mate is in storage.
							StorageIteratorType mate = findInStorage(_alignmentStorage, alignment->name());
							if(mate == _alignmentStorage.end()){
								//no mate found
								if(alignment->isProperPair()){
									//proper pair: add to storage and wait for mate
									addToContainer(_alignmentStorage, alignment, false);
								} else {
									//improper pair: add to blacklist and ready to write
									_blacklist.add(alignment->name());
									addToContainer(_alignmentStorage, alignment, true);
								}
							} else {
								if(alignment->readGroupId() != mate->alignment().readGroupId()){
									UERROR("Mates '", alignment->name(), "' are in different read groups!");
								}
								//mate found
								_handleMates(*alignment, mate);
							}
						}
					} else {
						//read is single end
						_handleSingle(*alignment);
					}
				}
			} else {
				//Did not pass QC: filter out
				//need to store in blacklist if it was paired
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
		_bamFile.printSummary(_outputName);
		_bamFile.close();
		_outBam.close();
	}
};

//-----------------------------------------
// TBamFilter
//-----------------------------------------
class TBamFilter final:public TGenomeParsedWithAlignmentStorage<TAlignmentStorage, TAlignmentStorageIterator>{
protected:
	void _openBamFileForWriting() override;
	void _handleMates(BAM::TAlignment & alignment, TAlignmentStorageIterator mate) override;
	void _handleSingle(BAM::TAlignment & alignment) override;
	bool _alignmentCanBeWrittenUnchanged() override;

public:
	TBamFilter() : TGenomeParsedWithAlignmentStorage() {}
	void run() {
		traverseBAM();
	}
};
} // namespace BamFilter
} // namespace GenomeTasks

#endif /* TBAMFILTER_H_ */
