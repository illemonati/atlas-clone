#include "TGenotypeMap.h"

//TGenotypeMap
TGenotypeMap::TGenotypeMap(){
    //create genotype map
    genotypeMap = new Genotype*[4];
    for(int i=0; i<4; ++i)
        genotypeMap[i] = new Genotype[4];

    //fill genotype map
    int geno = 0;
    for(int i=0; i<4; ++i){
        for(int j=i; j<4; ++j){
            genotypeMap[i][j] = static_cast<Genotype>(geno);
            genotypeMap[j][i] = genotypeMap[i][j];
            ++geno;
        }
    }

    //create and fill genotypeToBase
    numGenotypes = 10;
    genotypeToBase = new Base*[numGenotypes];
    for(int i=0; i<10; ++i){
        genotypeToBase[i] = new Base[2];
    }
    genotypeToBase[0][0] = A; genotypeToBase[0][1] = A;
    genotypeToBase[1][0] = A; genotypeToBase[1][1] = C;
    genotypeToBase[2][0] = A; genotypeToBase[2][1] = G;
    genotypeToBase[3][0] = A; genotypeToBase[3][1] = T;
    genotypeToBase[4][0] = C; genotypeToBase[4][1] = C;
    genotypeToBase[5][0] = C; genotypeToBase[5][1] = G;
    genotypeToBase[6][0] = C; genotypeToBase[6][1] = T;
    genotypeToBase[7][0] = G; genotypeToBase[7][1] = G;
    genotypeToBase[8][0] = G; genotypeToBase[8][1] = T;
    genotypeToBase[9][0] = T; genotypeToBase[9][1] = T;

    //create and fill context map
    numContexts = 25;
    numContextsNotN = 20;
    contextMap = new BaseContext*[5];

    //now fill regular context
    int context = 0;
    for(int i=0; i<5; ++i){
        contextMap[i] = new BaseContext[5];
        for(int j=0; j<4; ++j){
            contextMap[i][j] = static_cast<BaseContext>(context);
            ++context;
        }
    }

    //Now add those that should not occur, but sometimes do in bam files
    //Note that these should never occur in our data processing as they imply the base is N
    contextMap[0][4] = cAN;
    contextMap[1][4] = cCN;
    contextMap[2][4] = cGN;
    contextMap[3][4] = cTN;
    contextMap[4][4] = cNN;

    //fill base to char map
    baseToChar = new char[5];
    baseToChar[A] = 'A';
    baseToChar[C] = 'C';
    baseToChar[G] = 'G';
    baseToChar[T] = 'T';
    baseToChar[N] = 'N';

    //fill baseToFlippedBase map
    baseToFlippedBase = new Base[5];
    baseToFlippedBase[A] = T;
    baseToFlippedBase[C] = G;
    baseToFlippedBase[G] = C;
    baseToFlippedBase[T] = A;
    baseToFlippedBase[N] = N;

    //fill alleleicCombinationToBase
    alleleicCombinationToBase = new Base*[6];
    alleleicCombinationToBase[0] = new Base[2]; alleleicCombinationToBase[0][0] = A; alleleicCombinationToBase[0][1] = C;
    alleleicCombinationToBase[1] = new Base[2]; alleleicCombinationToBase[1][0] = A; alleleicCombinationToBase[1][1] = G;
    alleleicCombinationToBase[2] = new Base[2]; alleleicCombinationToBase[2][0] = A; alleleicCombinationToBase[2][1] = T;
    alleleicCombinationToBase[3] = new Base[2]; alleleicCombinationToBase[3][0] = C; alleleicCombinationToBase[3][1] = G;
    alleleicCombinationToBase[4] = new Base[2]; alleleicCombinationToBase[4][0] = C; alleleicCombinationToBase[4][1] = T;
    alleleicCombinationToBase[5] = new Base[2]; alleleicCombinationToBase[5][0] = G; alleleicCombinationToBase[5][1] = T;

    //fill alleleicCombinationToBaseChar
    alleleicCombinationToBaseChar = new char*[6];
    alleleicCombinationToBaseChar[0] = new char[2]; alleleicCombinationToBaseChar[0][0] = 'A'; alleleicCombinationToBaseChar[0][1] = 'C';
    alleleicCombinationToBaseChar[1] = new char[2]; alleleicCombinationToBaseChar[1][0] = 'A'; alleleicCombinationToBaseChar[1][1] = 'G';
    alleleicCombinationToBaseChar[2] = new char[2]; alleleicCombinationToBaseChar[2][0] = 'A'; alleleicCombinationToBaseChar[2][1] = 'T';
    alleleicCombinationToBaseChar[3] = new char[2]; alleleicCombinationToBaseChar[3][0] = 'C'; alleleicCombinationToBaseChar[3][1] = 'G';
    alleleicCombinationToBaseChar[4] = new char[2]; alleleicCombinationToBaseChar[4][0] = 'C'; alleleicCombinationToBaseChar[4][1] = 'T';
    alleleicCombinationToBaseChar[5] = new char[2]; alleleicCombinationToBaseChar[5][0] = 'G'; alleleicCombinationToBaseChar[5][1] = 'T';

    //fill alleleicCombinationToGenotypes (hom, het, hom2)
    alleleicCombinationToGenotypes = new Genotype*[6];
    alleleicCombinationToGenotypes[0] = new Genotype[3]; alleleicCombinationToGenotypes[0][0] = AA; alleleicCombinationToGenotypes[0][1] = AC; alleleicCombinationToGenotypes[0][2] = CC;
    alleleicCombinationToGenotypes[1] = new Genotype[3]; alleleicCombinationToGenotypes[1][0] = AA; alleleicCombinationToGenotypes[1][1] = AG; alleleicCombinationToGenotypes[1][2] = GG;
    alleleicCombinationToGenotypes[2] = new Genotype[3]; alleleicCombinationToGenotypes[2][0] = AA; alleleicCombinationToGenotypes[2][1] = AT; alleleicCombinationToGenotypes[2][2] = TT;
    alleleicCombinationToGenotypes[3] = new Genotype[3]; alleleicCombinationToGenotypes[3][0] = CC; alleleicCombinationToGenotypes[3][1] = CG; alleleicCombinationToGenotypes[3][2] = GG;
    alleleicCombinationToGenotypes[4] = new Genotype[3]; alleleicCombinationToGenotypes[4][0] = CC; alleleicCombinationToGenotypes[4][1] = CT; alleleicCombinationToGenotypes[4][2] = TT;
    alleleicCombinationToGenotypes[5] = new Genotype[3]; alleleicCombinationToGenotypes[5][0] = GG; alleleicCombinationToGenotypes[5][1] = GT; alleleicCombinationToGenotypes[5][2] = TT;
}

TGenotypeMap::~TGenotypeMap(){
    for(int i=0; i<4; ++i){
        delete[] genotypeMap[i];
    }
    for(int i=0; i<5; ++i){
        delete[] contextMap[i];
    }
    delete[] genotypeMap;
    for(int i=0; i<10; ++i)
        delete[] genotypeToBase[i];
    delete[] genotypeToBase;
    delete[] contextMap;
    delete[] baseToChar;
    delete[] baseToFlippedBase;

    for(int i=0; i<6; ++i){
        delete[] alleleicCombinationToBase[i];
        delete[] alleleicCombinationToBaseChar[i];
        delete[] alleleicCombinationToGenotypes[i];
    }
    delete[] alleleicCombinationToBase;
    delete[] alleleicCombinationToBaseChar;
    delete[] alleleicCombinationToGenotypes;
}

Base TGenotypeMap::getBase(const char base){
    if(base == 'A') return A;
    if(base == 'C') return C;
    if(base == 'G') return G;
    if(base == 'T') return T;
    if(base == 'a') return A;
    if(base == 'c') return C;
    if(base == 'g') return G;
    if(base == 't') return T;
    return N;
}

Base TGenotypeMap::getBaseOnlyCapitals(const char base){
    if(base == 'A') return A;
    if(base == 'C') return C;
    if(base == 'G') return G;
    if(base == 'T') return T;
    return N;
}

char TGenotypeMap::getBaseAsChar(Base base){
    if(base == A) return 'A';
    if(base == C) return 'C';
    if(base == G) return 'G';
    if(base == T) return 'T';
    return 'N';
}

char TGenotypeMap::getBaseAsChar(int base){
    if(base == A) return 'A';
    if(base == C) return 'C';
    if(base == G) return 'G';
    if(base == T) return 'T';
    return 'N';
}

Base TGenotypeMap::flipBase(char &base){
    if(base == 'A') return T;
    if(base == 'C') return G;
    if(base == 'G') return C;
    if(base == 'T') return A;
    if(base == 'a') return T;
    if(base == 'c') return G;
    if(base == 'g') return C;
    if(base == 't') return A;
    return N;
}

Genotype TGenotypeMap::getGenotype(Base first, Base second){
    if(first < second) return genotypeMap[first][second];
    else return genotypeMap[second][first];
}

Genotype TGenotypeMap::getGenotype(int first, int second){
    if(first < second) return genotypeMap[first][second];
    else return genotypeMap[second][first];
}

Genotype TGenotypeMap::getGenotype(char first, char second){
    Base Bfirst = getBase(first);
    Base Bsecond = getBase(second);
    if(Bfirst < Bsecond) return genotypeMap[Bfirst][Bsecond];
    else return genotypeMap[Bsecond][Bfirst];
}

std::string TGenotypeMap::getGenotypeString(int num){
    if(num==0) return "AA";
    if(num==1) return "AC";
    if(num==2) return "AG";
    if(num==3) return "AT";
    if(num==4) return "CC";
    if(num==5) return "CG";
    if(num==6) return "CT";
    if(num==7) return "GG";
    if(num==8) return "GT";
    if(num==9) return "TT";
    throw "GenotypeMap: Unknown genotype with number " + toString(num) + "!";
}

std::string TGenotypeMap::getGenotypeString(Genotype geno){
    if(geno==0) return "AA";
    if(geno==1) return "AC";
    if(geno==2) return "AG";
    if(geno==3) return "AT";
    if(geno==4) return "CC";
    if(geno==5) return "CG";
    if(geno==6) return "CT";
    if(geno==7) return "GG";
    if(geno==8) return "GT";
    if(geno==9) return "TT";
    throw "GenotypeMap: Unknown genotype with number " + toString((int) geno) + "!";
}

std::string TGenotypeMap::getGenotypeStringKnownAlleles(int num, char ref, char alt){
    std::string genotype = "";
    if(num == 0){
        genotype += ref;
        genotype += ref;
        return genotype;
    }
    else if(num == 2){
        genotype += alt;
        genotype += alt;
        return genotype;
    }
    else if(num == 1){
        Base refBase = getBase(ref);
        Base altBase = getBase(alt);
        if(refBase > altBase){
            genotype += alt;
            genotype += ref;
            return genotype;
        }
        else{
            genotype += ref;
            genotype += alt;
            return genotype;
        }
    }
    throw "GenotypeMap: Unknown genotype with number " + toString(num) + "!";
}

std::pair<Base, Base> TGenotypeMap::getBasesOfGenotype(int num){
    if(num==0) return std::pair<Base,Base>(A,A);
    if(num==1) return std::pair<Base,Base>(A,C);
    if(num==2) return std::pair<Base,Base>(A,G);
    if(num==3) return std::pair<Base,Base>(A,T);
    if(num==4) return std::pair<Base,Base>(C,C);
    if(num==5) return std::pair<Base,Base>(C,G);
    if(num==6) return std::pair<Base,Base>(C,T);
    if(num==7) return std::pair<Base,Base>(G,G);
    if(num==8) return std::pair<Base,Base>(G,T);
    if(num==9) return std::pair<Base,Base>(T,T);
    throw "GenotypeMap: Unknown genotype with number " + toString(num) + "!";
}

int TGenotypeMap::getNumContext(){
    return 20;
}


//TBaseContext
BaseContext TGenotypeMap::getContext(Base first, Base second){
    if(second == N) throw "Context not defined with second base = N!";
    return contextMap[first][second];
}

BaseContext TGenotypeMap::getContext(int first, int second){
    if(second > 3) throw "Context not defined with second base = N!";
    return contextMap[first][second];
}

BaseContext TGenotypeMap::getContext(char first, char second){
    return getContext(getBase(first), getBase(second));
}

BaseContext TGenotypeMap::getContextReverseRead(char first, char second){
    return getContext(flipBase(first), flipBase(second));
}

std::string TGenotypeMap::getContextString(int num){
    if(num==0) return "AA";
    if(num==1) return "AC";
    if(num==2) return "AG";
    if(num==3) return "AT";
    if(num==4) return "CA";
    if(num==5) return "CC";
    if(num==6) return "CG";
    if(num==7) return "CT";
    if(num==8) return "GA";
    if(num==9) return "GC";
    if(num==10) return "GG";
    if(num==11) return "GT";
    if(num==12) return "TA";
    if(num==13) return "TC";
    if(num==14) return "TG";
    if(num==15) return "TT";
    if(num==16) return "-A";
    if(num==17) return "-C";
    if(num==18) return "-G";
    if(num==19) return "-T";
    throw "GenotypeMap: Unknown text with number " + toString(num) + "!";
}

TBaseFrequencies::TBaseFrequencies(){
    for(int i = 0; i < 4; ++i) freq[i] = 0.0;
    wasNormalized = false;
}

void TBaseFrequencies::add(Base B, double &weight){
    freq[B] += weight;
}

void TBaseFrequencies::addNoRef(Base B, double weight){
    freq[B] += weight;
}

void TBaseFrequencies::normalize(){
    if(!wasNormalized){
        double sum = 0.0;
        for(int i = 0; i < 4; ++i) sum += freq[i];
        sum += 4.0;
        for(int i = 0; i < 4; ++i) freq[i] = (freq[i] + 1.0) / sum;
        wasNormalized = true;
    }
}

void TBaseFrequencies::setEqualBaseFreq(){
    for(int i = 0; i < 4; ++i) freq[i] = 0.25;
}

void TBaseFrequencies::clear(){
    for(int i = 0; i < 4; ++i) freq[i] = 0.0;
    wasNormalized = false;
}

void TBaseFrequencies::print(){
    std::cout << "freq(A) = " << freq[0] << ", freq(C) = " << freq[1] << ", freq(G) = " << freq[2] << ", freq(T) = " << freq[3] << std::endl;
}

double &TBaseFrequencies::operator[](int pos){
    return freq[pos];
}
