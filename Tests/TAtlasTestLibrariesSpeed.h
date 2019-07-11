#ifndef TATLASTESTLIBRARIESSPEED_H
#define TATLASTESTLIBRARIESSPEED_H

#include <vector>
#include <iterator>
#include <string>
#include <chrono>

#include "IOAbstractClasses/BamReader.h"
#include "IOAbstractClasses/BamAlignment.h"
#include "IOAbstractClasses/Fasta.h"
#include "IOTool.h"

#include "TAtlasTest.h"
#include "TParameters.h"
#include "TLog.h"


#include <iomanip>


class TAtlasTestLibrarieSpeed : public TAtlasTest
{
private:
    std::string bamFile;
    std::string fastaFile;
    std::vector<BamReader*> bamReaders;
    std::vector<Fasta*> fastaReaders;
    std::vector<std::string> librariesName;
    TLog* logfile;

    void testSpeedBamReading();
    void testSpeedFastaReading();
    std::string convertNanoSecondToReadableString(long ms);

public:
    TAtlasTestLibrarieSpeed(TParameters & params, TLog* logfile);
    bool run();
};

#endif // TATLASTESTLIBRARIESSPEED_H
