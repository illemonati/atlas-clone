#make file for atlas

SRC = $(wildcard *.cpp) $(wildcard *.C) $(wildcard Tests/*.cpp) $(wildcard Simulations/*.cpp) $(wildcard Vcf/*.cpp) $(wildcard Recalibration/*.cpp) $(wildcard PopulationTools/*.cpp) $(wildcard IOTools/*.cpp)  $(wildcard IOTools/BamToolsAdapters/*.cpp) $(wildcard IOTools/Compression/*.cpp) $(wildcard Utils/*.cpp) $(wildcard Alignment/*.cpp) $(wildcard Bed/*.cpp) $(wildcard Site/*.cpp) $(wildcard Task/*.cpp) $(wildcard Theta/*.cpp) $(wildcard bamtools/api/*.cpp) $(wildcard bamtools/api/algorithms/*.cpp) $(wildcard bamtools/api/internal/bam/*.cpp) $(wildcard bamtools/api/internal/index/*.cpp) $(wildcard bamtools/api/internal/io/*.cpp) $(wildcard bamtools/api/internal/sam/*.cpp) $(wildcard bamtools/api/internal/utils/*.cpp) $(wildcard bamtools/utils/*.cpp)

GIT_HEADER = gitversion.cpp
OBJ = $(SRC:%.cpp=%.o)
BIN = atlas

INCLUDE = -I. -Ibamtools -ITests -ISimulations -IVcf -IRecalibration -IPopulationTools -IIOTools -IIOTools/Compression -IUtils -IAlignment -IBed -ISite -ITask -ITheta

.PHONY : all
all : $(BIN)

#If armadillo not available. Default true 
ifeq ($(ARM),)
BINFLAG = -lz -larmadillo 
OBJFLAG = -std=c++1y
else
BINFLAG = -lz -lblas -llapack
OBJFLAG = -Iarmadillo-8.400.0/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -std=c++1y
endif

#If SeqLib not available. Default true 
ifeq ($(SEQLIB),)
DEFINE = -DSEQLIB
LIB = SeqLib/lib/libseqlib.a SeqLib/lib/libbwa.a SeqLib/lib/libhts.a SeqLib/lib/libfml.a
BINFLAG += -lpthread -llzma -lbz2 
INCLUDE += -ISeqLib/htslib -ISeqLib
SRC += $(wildcard IOTools/SeqLibAdapters/*.cpp)
endif

#If we want to debug. Default false
ifeq ($(DEBUG),true)
DEFINE += -DDEBUG
DEBUGFLAG = -g -lprofiler
endif

#If LZ4 not available. Default true 
ifeq ($(LZ4),)
DEFINE += -DLZ4
BINFLAG += -llz4
endif

#If ZSTD not available. Default true
ifeq ($(ZSTD),)
DEFINE += -DZSTD
BINFLAG += -lzstd
endif

$(BIN): $(GIT_HEADER) $(OBJ)
	$(CXX) -rdynamic -O3 $(DEBUGFLAG) -o $(BIN) $(OBJ) $(LIB) $(BINFLAG)


$(GIT_HEADER): .git/HEAD .git/COMMIT_EDITMSG
	echo "#include \"gitversion.h\"" > $@
	echo "std::string getGitVersion(){" >> $@
	echo "return \"$(shell git rev-parse HEAD)\";" >> $@
	echo "}" >> $@

.git/COMMIT_EDITMSG :
	touch $@


%.o: %.cpp
	$(CXX) -O3 $(DEBUGFLAG) -c $(DEFINE) $(INCLUDE) $(OBJFLAG) $< -o $@


.PHONY : clean
clean:
	rm -rf $(BIN) $(OBJ)
