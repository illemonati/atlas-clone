/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TMain.h"
//#include <gperftools/profiler.h>

//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
//    ProfilerStart("majorminor_2_big_sample.prof");
	TMain main("Atlas", "0.9", "https://bitbucket.org/wegmannlab/atlas");
    bool r = main.run(argc, argv);
//    ProfilerStop();
    return r;
};

