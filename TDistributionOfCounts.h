/*
 * TDepthCounts.h
 *
 *  Created on: Mar 18, 2019
 *      Author: phaentu
 */

#ifndef TDISTRIBUTIONOFCOUNTS_H_
#define TDISTRIBUTIONOFCOUNTS_H_

#include "TFile.h"

class TDistributionOfCounts{
private:
	int maxValue;
	int numBins;
	int lastBin;
	std::string whatIsCounted;

	unsigned int* counts;

public:
	TDistributionOfCounts(int MaxValue, std::string nameOfCounts){
		whatIsCounted = nameOfCounts;
		maxValue = MaxValue;
		numBins = maxValue + 2; //starts at zero and includes bin to lump > maxValue
		lastBin = maxValue + 1;
		counts = new unsigned int[numBins];
		for(int i=0; i<numBins; ++i)
			counts[i] = 0;
	};

	~TDistributionOfCounts(){
		delete[] counts;
	};

	void add(int depth){
		if(depth > maxValue)
			++counts[lastBin];
		else
			++counts[depth];
	};

	void add(int depth, int number){
		if(depth > maxValue)
			counts[lastBin] += number;
		else
			counts[depth] += number;
	};

	unsigned int totNumSites(){
		unsigned int sum = 0;
		for(int i=0; i<numBins; ++i)
			sum += counts[i];
		return sum;
	};

	void writeCounts(std::string filename){
		//open file
		TOutputFile out(filename);
		out.writeHeader({whatIsCounted, "counts"});

		//write bins
		for(int i=0; i<=maxValue; i++)
			out << i << counts[i] << std::endl;

		//write last bin
		out << ">" + toString(maxValue) << counts[lastBin] << std::endl;
	};

	void writeNormalizedCounts(std::string filename){
		//open file
		TOutputFile out(filename);
		out.writeHeader({whatIsCounted, "density"});

		//get sum
		double sum = totNumSites();

		//write bins
		for(int i=0; i<=maxValue; i++)
			out << i << counts[i] / sum << std::endl;

		//write last bin
		out << ">" + toString(maxValue) << counts[lastBin] / sum << std::endl;
	};

	void writeCumulativeCounts(std::string filename){
		//open file
		TOutputFile out(filename);
		out.writeHeader({whatIsCounted, "cumulativeCounts"});

		//write bins
		unsigned cumul = 0;
		for(int i=0; i<=maxValue; i++){
			cumul += counts[i];
			out << i << cumul << std::endl;
		}

		//write last bin
		out << ">" + toString(maxValue) << cumul + counts[lastBin] << std::endl;
	};

	void writeNormalizedCumulativeCounts(std::string filename){
		//open file
		TOutputFile out(filename);
		out.writeHeader({whatIsCounted, "cumulativeDensity"});

		//get sum
		double sum = totNumSites();

		//write bins
		double cumul = 0.0;
		for(int i=0; i<=maxValue; i++){
			cumul += counts[i] / sum;
			out << i << cumul << std::endl;
		}

		//write last bin
		out << ">" + toString(maxValue) << 1.0 << std::endl;
	};

	void writeQuantiles(std::string filename){
		//open file
		TOutputFile out(filename);
		out.writeHeader({"quantile", whatIsCounted});

		//normalized cumulative distribution and quantiles
		float quantiles[16] = {0.0, 0.001, 0.005, 0.01, 0.025, 0.05, 0.2, 0.5, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 1.0};
		int numQuantiles = 16;

		//get sum
		double sum = totNumSites();

		//find quantiles
		double cumul = 0.0;
		int nextQuantile = 0;

		for(int i=0; i<=maxValue; i++){
			cumul += counts[i] / sum;

			while(cumul >= quantiles[nextQuantile] && nextQuantile < numQuantiles){
				out << quantiles[nextQuantile] << i << std::endl;
				++nextQuantile;
			}
		}

		//write >max for remaining
		for(int p=nextQuantile; p<numQuantiles; ++p){
			out << quantiles[p] << ">" + toString(maxValue) << std::endl;
		}
	};
};


#endif /* TDISTRIBUTIONOFCOUNTS_H_ */
