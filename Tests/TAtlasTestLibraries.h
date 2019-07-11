#ifndef TATLASTESTLIBRARIES_H
#define TATLASTESTLIBRARIES_H

#include "IOTool.h"
#include "IOAbstractClasses/Fasta.h"
#include "IOAbstractClasses/BamReader.h"
#include "IOAbstractClasses/BamAlignment.h"
#include "IOAbstractClasses/SamHeader.h"
#include "IOAbstractClasses/RefVector.h"

#include "TAtlasTest.h"
#include "TParameters.h"
#include "TLog.h"
#include <unordered_set>
#include <string>


class TAtlasTestLibrarieCompareFunctions : public TAtlasTest
{
private:
    std::string bamFile;
    std::string fastaFile;
    std::string lib1;
    std::string lib2;
    TLog* logfile;

    //See if the two librairies have the same behaviour when reading a bam file
    bool compareBamReading();
    //See if the two librairies have the same behaviour when reading a fasta file
    bool compareFastaReading();
    //Check if the result is the same for string
    void checkResultMatchStr(std::string function, std::string resLib1, std::string resLib2);
    //Check if the result is the same for integer
    void checkResultMatchInt(std::string function, int resLib1, int resLib2);
    //Write function mismatching
    void writeMismatch(std::string function, std::string resLib1, std::string resLib2);

public:
    TAtlasTestLibrarieCompareFunctions(TParameters & params, TLog* logfile);
    bool run();

};

#endif // TATLASTESTLIBRARIES_H
