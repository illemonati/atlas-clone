#include "TReadGroup.h"
//---------------------------------------------------------------
//TReadGroup
//---------------------------------------------------------------
namespace BAM {

std::string TReadGroup::compileSamHeader() const{
	std::string header = "@RG\tID:" + name_ID;

	if(!barcodeSequence_BC.empty()){
		header += "\tBC:" + barcodeSequence_BC;
	}

	if(!sequencingCenter_CN.empty()){
		header += "\tCN:" + sequencingCenter_CN;
	}

	if(!description_DS.empty()){
		header += "\tDS:" + description_DS;
	}

	if(!productionDate_DT.empty()){
		header += "\tDT:" + productionDate_DT;
	}

	if(!flowOrder_FO.empty()){
		header += "\tFO:" + flowOrder_FO;
	}

	if(!keySequence_KS.empty()){
		header += "\tKS:" + keySequence_KS;
	}

	if(!library_LB.empty()){
		header += "\tLB:" + library_LB;
	}

	if(!program_PG.empty()){
		header += "\tPG:" + program_PG;
	}

	if(!predictedInsertSize_PI.empty()){
		header += "\tPI:" + predictedInsertSize_PI;
	}

	if(!sequencingTechnology_PL.empty()){
		header += "\tPL:" + sequencingTechnology_PL;
	}

	if(!platformModel_PM.empty()){
		header += "\tPM:" + platformModel_PM;
	}

	if(!platformUnit_PU.empty()){
		header += "\tPU:" + platformUnit_PU;
	}

	if(!sample_SM.empty()){
		header += "\tSM:" + sample_SM;
	}

	return header;
}

}
