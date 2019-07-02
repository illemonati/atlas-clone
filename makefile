#make file for atlas

SRC = $(wildcard *.cpp) $(wildcard *.C) $(wildcard Tests/*.cpp) $(wildcard Simulations/*.cpp) $(wildcard Vcf/*.cpp) $(wildcard Recalibration/*.cpp) $(wildcard PopulationTools/*.cpp) $(wildcard IOTools/*.cpp) $(wildcard IOTools/SeqLibAdapters/*.cpp) $(wildcard IOTools/BamToolsAdapters/*.cpp) $(wildcard Utils/*.cpp) $(wildcard Alignment/*.cpp) $(wildcard Bed/*.cpp) $(wildcard Site/*.cpp) $(wildcard Task/*.cpp) $(wildcard Theta/*.cpp) $(wildcard bamtools/api/*.cpp) $(wildcard bamtools/api/algorithms/*.cpp) $(wildcard bamtools/api/internal/bam/*.cpp) $(wildcard bamtools/api/internal/index/*.cpp) $(wildcard bamtools/api/internal/io/*.cpp) $(wildcard bamtools/api/internal/sam/*.cpp) $(wildcard bamtools/api/internal/utils/*.cpp) $(wildcard bamtools/utils/*.cpp)

GIT_HEADER = gitversion.cpp
OBJ = $(SRC:%.cpp=%.o)
BIN = atlas
LIB = SeqLib/lib/libseqlib.a SeqLib/lib/libbwa.a SeqLib/lib/libhts.a SeqLib/lib/libfml.a
INCLUDE = -I. -Ibamtools -ITests -ISimulations -IVcf -IRecalibration -IPopulationTools -IIOTools -IUtils -IAlignment -IBed -ISite -ITask -ITheta -ISeqLib/htslib -ISeqLib

.PHONY : all
all : $(BIN)

ifeq ($(ARM),)
BINFLAG = -lz -larmadillo -lpthread -llzma -lbz2
OBJFLAG = -std=c++1y
else
BINFLAG = -lz -lblas -llapack -lpthread -llzma -lbz2
OBJFLAG = -Iarmadillo-8.400.0/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -std=c++1y
endif

$(BIN): $(GIT_HEADER) $(OBJ)
	$(CXX) -O3 -g -lprofiler -o $(BIN) $(OBJ) $(LIB) $(BINFLAG)


$(GIT_HEADER): .git/HEAD .git/COMMIT_EDITMSG
	echo "#include \"gitversion.h\"" > $@
	echo "std::string getGitVersion(){" >> $@
	echo "return \"$(shell git rev-parse HEAD)\";" >> $@
	echo "}" >> $@

.git/COMMIT_EDITMSG :
	touch $@


%.o: %.cpp
	$(CXX) -O3 -g -lprofiler -c $(INCLUDE) $(OBJFLAG) $< -o $@


.PHONY : clean
clean:
	rm -rf $(BIN) $(OBJ)
