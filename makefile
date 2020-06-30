#make file for atlas

SRC = $(wildcard *.cpp) $(wildcard *.C) $(wildcard Tests/*.cpp) $(wildcard commonutilities/*.cpp) $(wildcard commonutilities/IntegrationTests/*.cpp) $(wildcard Simulations/*.cpp) $(wildcard Vcf/*.cpp) $(wildcard Recalibration/*.cpp) $(wildcard PopulationTools/*.cpp) $(wildcard GLF/*.cpp) $(wildcard bamtools/api/*.cpp) $(wildcard bamtools/api/algorithms/*.cpp) $(wildcard bamtools/api/internal/bam/*.cpp) $(wildcard bamtools/api/internal/index/*.cpp) $(wildcard bamtools/api/internal/io/*.cpp) $(wildcard bamtools/api/internal/sam/*.cpp) $(wildcard bamtools/api/internal/utils/*.cpp) $(wildcard bamtools/utils/*.cpp)
GIT_HEADER = commonutilities/gitversion.cpp

# check if commonutilities/gitversion.cpp already exists (if it doesn't, we need to add it to SRC)
ifeq (,$(findstring $(GIT_HEADER),$(SRC)))
    SRC += $(GIT_HEADER)
endif

OBJ = $(SRC:%.cpp=%.o)
BIN = atlas

.PHONY : all
all : $(BIN)

# create gitversion.cpp file
$(GIT_HEADER): .git/HEAD .git/COMMIT_EDITMSG .git/index
	touch $@
	echo "#include \"gitversion.h\"" > $@
	echo "std::string getGitVersion(){" >> $@
	echo "return \"$(shell git rev-parse HEAD)\";" >> $@
	echo "}" >> $@


ifeq ($(ARM),)
BINFLAG = -lz -larmadillo
OBJFLAG = -std=c++1y
else
BINFLAG = -lz -lblas -llapack
OBJFLAG = -Iarmadillo-8.400.0/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -std=c++1y
endif

$(BIN): $(GIT_HEADER) $(OBJ)
	$(CXX) -O3 -o $(BIN) $(OBJ) $(BINFLAG)

.git/COMMIT_EDITMSG :
	touch $@


%.o: %.cpp
	$(CXX) -O3 -c -I. -Ibamtools -Icommonutilities -IGLF -ITests -ISimulations -IVcf -IRecalibration -IPopulationTools $(OBJFLAG)  $< -o $@


.PHONY : clean
clean:
	rm -rf $(BIN) $(OBJ)
