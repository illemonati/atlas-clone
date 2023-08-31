/*
 * TPileup.cpp
 *
 *  Created on: Jun 4, 2020
 *      Author: wegmannd
 */

#include "TPileup.h"

#include <stdint.h>
#include <algorithm>
#include <ostream>
#include <vector>

#include "genometools/GenotypeTypes.h"
#include "TGenotypeData.h"
#include "TGenotypeLikelihoodCalculator.h"
#include "coretools/Main/TLog.h"
#include "coretools/Main/TParameters.h"
#include "TSite.h"
#include "TWindow.h"

namespace GenomeTasks{

using coretools::instances::logfile;
using coretools::instances::parameters;

namespace impl {

void addNames(std::vector<std::string> &vec) {
	using GT = genometools::Genotype;
	for (auto g = GT::min; g < GT::max; ++g) vec.push_back("P(D|" + genometools::toString(g) + ")");
};

void write(const GenotypeLikelihoods::TGenotypeLikelihoods& lhs, coretools::TOutputFile & out) {
	for (const auto l: lhs)
		out << l;
};

bool parseField(std::set<std::string> &fields, const std::string &tag, const std::string &explanation) {
	if (fields.find(tag) != fields.end()) {
		logfile().list(explanation + " (" + tag + ")");
		fields.erase(fields.find(tag));
		return true;
	}
	return false;
}

} // namespace

//---------------------------------
// TPileup
//---------------------------------
TPileup::TPileup():TGenome_windows(){
	if(!parameters().parameterExists("onlySummaries")){
		//open output file
		const std::string filename = _outputName + "_pileup.txt.gz";
		logfile().list("Writing pileup to file '" + filename + "'.");
		_out.open(filename);

		//parse output fields
		logfile().startIndent("Will print the following pileup fields (parameter 'fields'):");
		std::set<std::string> fields;
		parameters().fillParameterIntoContainerWithDefault("fields", fields, ',', {"depth", "bases", "qualities", "alleles", "mates", "strands", "likelihoods"});

		_printSettings.set<Depth>(impl::parseField(fields, "depth", "Sequencing depth"));
		_printSettings.set<Bases>(impl::parseField(fields, "bases", "Pileup bases"));
		_printSettings.set<Qualities>(impl::parseField(fields, "qualities", "Pileup qualities"));
		_printSettings.set<Alleles>(impl::parseField(fields, "alleles", "Allele counts"));
		_printSettings.set<Mates>(impl::parseField(fields, "mates", "Mate information"));
		_printSettings.set<Strand>(impl::parseField(fields, "strands", "Strand information"));
		_printSettings.set<Likelihoods>(impl::parseField(fields, "likelihoods", "Genotype likelihoods"));
		logfile().endIndent();

		//check if unknown fields were given
		if(!fields.empty()){
			if(fields.size() == 1){
				UERROR("Unknown field '", *fields.begin(), "'! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates' and 'strands'.");
			} else {
				std::string f;
				for(auto i : fields){
					f += '\'' + i + "', ";
				}
				UERROR("Unknown fields: ", f.substr(0, f.size() - 2), "! Valid fields are 'depth', 'bases', 'qualities', 'alleles', 'mates',  'strands' and 'likelihoods'.");
			}
		}

		//compile header
		std::vector<std::string> header = {"chromosome", "position"};
		if(_reference){ header.push_back("reference"); }
		if(_printSettings.get<Depth>()){
			header.push_back("depth");
			if(_reference){
				header.push_back("refDepth");
			}
		}
		if(_printSettings.get<Bases>()){ header.push_back("bases"); }
		if(_printSettings.get<Qualities>()){ header.push_back("qualities"); }
		if(_printSettings.get<Alleles>()){
			header.push_back("numA");
			header.push_back("numC");
			header.push_back("numG");
			header.push_back("numT");
			if(_reference){
				header.push_back("numRef");
				header.push_back("numNonRef");
			}
			header.push_back("numAlleles");
		}
		if(_printSettings.get<Mates>()){
			header.push_back("numFirstMate");
			header.push_back("numSecondMate");
		}
		if(_printSettings.get<Strand>()){
			header.push_back("numForwardStrand");
			header.push_back("numReverseStrand");
		}
		if(_printSettings.get<Likelihoods>()){
			impl::addNames(header);
		}

		_out.writeHeader(header);

		//print all sites, also those without data?
		if(parameters().parameterExists("printAll")){
			_printSettings.set<OnlySitesWithData>(false);
			logfile().list("Will print all sites that pass filters, including those without data. (parameter 'printAll')");
		} else {
			_printSettings.set<OnlySitesWithData>(true);
			logfile().list("Will print only sites with data. (use 'printAll' to print all)");
		}
	}

	if (parameters().parameterExists("histograms") || parameters().parameterExists("onlySummaries")){
		logfile().startIndent("Will print the following histograms (parameter 'histograms'):");
		std::set<std::string> histograms;
		parameters().fillParameterIntoContainerWithDefault("histograms", histograms, ',', {"depth", "qualities", "contexts"});

		if (histograms.empty()){
			histograms.emplace("depth");
			histograms.emplace("qualities");
			histograms.emplace("contexts");
		}

		_histSettings.set<Depths>(impl::parseField(histograms, "depth", "Sequencing depth"));
		_histSettings.set<Quality>(impl::parseField(histograms, "qualities", "Base qualities"));
		_histSettings.set<Contexts>(impl::parseField(histograms, "contexts", "Base contexts"));
		_histSettings.set<AllelicDepth>(impl::parseField(histograms, "allelicDepth", "Allelic depth"));
		logfile().endIndent();

		//check if unknown fields were given
		if(!histograms.empty()){
			if(histograms.size() == 1){
				UERROR("Unknown histogram '", *histograms.begin(), "'! Valid histograms are 'depth', 'qualities', 'allelicDepth' and 'contexts'.");
			} else {
				std::string f;
				for(auto i : histograms){
					f += '\'' + i + "', ";
				}
				UERROR("Unknown histograms: ", f.substr(0, f.size() - 2), "! Valid histograms are 'depth', 'qualities', 'allelicDepth' and 'contexts'.");
			}
		}
		if(_histSettings.get<Depths>()){
			_outDepthHistogram.open(_outputName + "_depthPerWindow.txt.gz", {"window", "depth"});
			_outDepthPerChromosome.open(_outputName + "_depthPerChromosome.txt.gz", {"chromosome", "depth"});
		}

		if(_histSettings.get<AllelicDepth>()){
			logfile().list("Will assemble allelic depth up to a max depth of " + coretools::str::toString(_readUpToDepth) + ". (parameter 'readUpToDepth')");
			if(_readUpToDepth > 100){
				logfile().warning("Allocating count table for a max depth of " + coretools::str::toString(_readUpToDepth) + " uses a lot of memory! Use argument readUpToDepth to limit.");
			}

			_counts.resize(_readUpToDepth);

			if(parameters().parameterExists("includeZero")){
				_writeEmpty = true;
				logfile().list("Will write full table, including cells with zero counts. (parameter 'includeZero')");
			} else {
				_writeEmpty = false;
				logfile().list("Will only print cells with non-zero counts. (use 'includeZero' to print all cells)");
			}
		}
	} else {
		logfile().list("Will not output histograms (use 'histograms' to do so).");
	}
};


void TPileup::_handleWindow(){
	using genometools::Base;
	logfile().list("Writing pileup ...");

	if(_histSettings.get<AllelicDepth>()){
		logfile().list("Adding sites to allelic depth table ...");
	}
	for (size_t pos = 0; pos < _window.size(); ++pos) {
		const auto &site = _window[pos];
		if(!parameters().parameterExists("onlySummaries")){
			if (_printSettings.get<OnlySitesWithData>() && site.empty()) continue;
			_out << _window.chrName();
			_out << _window.positionOnChr(pos) + 1; // positions are zero-based internally

			if(_reference){
				_out << site.refBase;
			}
			if(_printSettings.get<Depth>()){
				_out << site.depth();
				if(_reference){
					_out << site.refDepth();
				}
			}
			if(_printSettings.get<Bases>()){
				_out << site.getBases();
			}
			if(_printSettings.get<Qualities>()){
				_out << site.getQualities();
			}
			if(_printSettings.get<Alleles>()){
				const auto alleleCounts = site.countAlleles();
				_out << alleleCounts[Base::A] << alleleCounts[Base::C] << alleleCounts[Base::G] << alleleCounts[Base::T];
				if(_reference){
					_out << alleleCounts[site.refBase] << alleleCounts.size() - alleleCounts[site.refBase];
				}
				_out << (int) coretools::numNonZero(alleleCounts);
			}
			if(_printSettings.get<Mates>()){
				const auto  mateCounts = site.countMates();
				_out << mateCounts[BAM::Mate::first] << mateCounts[BAM::Mate::second];
			}
			if(_printSettings.get<Strand>()){
				const auto  strandCounts = site.countFwdRev();
				_out << strandCounts[0] << strandCounts[1];
			}
			if(_printSettings.get<Likelihoods>()){
				_genoLik = _genotypeLikelihoodCalculator.calculateGenotypeLikelihoods(site);
				impl::write(_genoLik, _out);
			}
		_out.endln();
		}
		//write histograms
		if(_histSettings.get<Depths>()){
			_depthPerSite.add(site.depth());
			_depthPerSitePerChromosome.add(site.depth());
		}

		if(_histSettings.get<Quality>()){
			for (auto &b: site){
				if(b.base != genometools::Base::N){
					_qualDist.add(b.readGroupID, b.recalibratedQualityAsPhredInt.get());
				}
			}
		}

		if(_histSettings.get<Contexts>()){
			for(auto& b : site){
				_contextDist.add(b.recalibratedQualityAsPhredInt.get(), coretools::index(b.context()));
			}
		}

		if(_histSettings.get<AllelicDepth>()){
			const auto alleleCounts = site.countAlleles();
			_counts.addSite(alleleCounts);
		}
	}

	//write depth per window
	//also write depth per chromosome if this window is the last window of a chromosome
	if(_histSettings.get<Depths>()){
		logfile().list("Writing sequencing depth estimates to file ...");
		_outDepthHistogram.writeNoDelim(_window.chrName(), ':', _window.from().position() + 1, '-', _window.to().position()).writeDelim();
		_outDepthHistogram.writeln(_window.depth());
		logfile().done();
		if(_bamFile.chrChanged()){
                        _outDepthPerChromosome.writeln(_window.chrName(),_depthPerSitePerChromosome.mean());
                        _depthPerSitePerChromosome.clear();
                }
	}
};

void TPileup::run(){
	_traverseBAMWindows();

	if (_histSettings.get<Depths>()){
		//write distribution
		logfile().list("Writing depth per site distribution to file '", _outputName, "_depthPerSiteHistogram.txt.gz'"); 
		logfile().list("Writing average depth per window to file '", _outputName, "_depthPerWindow.txt.gz'");
		logfile().list("Writing average depth per chromosome to file '", _outputName, "_depthPerChromosome.txt.gz'");
		_depthPerSite.write(_outputName + "_depthPerSiteHistogram.txt.gz", "depth");
		_outDepthPerChromosome.writeln(_window.chrName(),_depthPerSitePerChromosome.mean());
	}

	if (_histSettings.get<Quality>()){
		//print distribution
		std::string outputFileName = _outputName + "_qualHistogram.txt.gz";
		logfile().list("Writing quality distribution to '" + outputFileName + "'.");
		coretools::TOutputFile out(outputFileName, {"readGroup", "quality", "counts"});

		//get read group names
		std::vector<std::string> readGroupNames;
		_bamFile.readGroups().fillVectorWithNames(readGroupNames);
		//write combined
		_qualDist.writeCombined(out, "allReadGroups");
		_qualDist.write(out, readGroupNames);	
	}

	if (_histSettings.get<Contexts>()){
		//write counts
		std::string outputFileName = _outputName + "_contextInformation.txt.gz";
		logfile().list("Writing context information to file '" + outputFileName + "'.");

		std::vector<std::string> contextLabels;

		for(genometools::BaseContext c = genometools::BaseContext::min; c < genometools::BaseContext::max; ++c){
			contextLabels.push_back(genometools::toString(c));
		}

		_contextDist.writeAsMatrix(outputFileName, "quality", contextLabels);
	}

	if (_histSettings.get<AllelicDepth>()){
		//write to file
		std::string outputFileName = _outputName + "_allelicDepth.txt.gz";
		logfile().list("Writing allelic depth table to '" + outputFileName + "' ...");
		_counts.write(outputFileName, _writeEmpty);
		logfile().done();
	}
};

}; // end namespace
