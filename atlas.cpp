/*
 * estimHet.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: wegmannd
 */

#include "TMain.h"
#ifdef DEBUG
#include <gperftools/profiler.h>
#endif
//---------------------------------------------------------------------------
//Main function
//---------------------------------------------------------------------------
int main(int argc, char* argv[]){
#ifdef DEBUG
    ProfilerStart("majorminor_2_big_sample.prof");
#endif

	TMain main("Atlas", "0.9", "https://bitbucket.org/wegmannlab/atlas");
    bool r = main.run(argc, argv);

#ifdef DEBUG
    ProfilerStop();
#endif

    return r;
};
