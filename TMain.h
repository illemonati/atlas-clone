/*
 * TMAin.h
 *
 *  Created on: Apr 1, 2019
 *      Author: phaentu
 */

#ifndef TMAIN_H_
#define TMAIN_H_

#include "TTaskList.h"
#include "gitversion.h"
#include "Tests/TTesting.h"
#include "IOTools/IOTool.h"

//---------------------------------------------------------------------------
// TMain
//---------------------------------------------------------------------------
class TMain{
private:
	struct timeval start;
	std::string applicationName;
	std::string version;
	std::string web;
	std::string title;
	TLog logfile;

    void fillTitleString();

public:
    TMain(std::string ApplicationName, std::string Version, std::string Web);
    int run(int argc, char* argv[]);
};


#endif /* TMAIN_H_ */
