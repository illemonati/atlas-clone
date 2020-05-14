#make file for atlas

SRC = $(wildcard *.cpp) $(wildcard *.C) $(wildcard Tests/*.cpp) $(wildcard commonutilities/*.cpp) $(wildcard commonutilities/IntegrationTests/*.cpp) $(wildcard Simulations/*.cpp) $(wildcard Vcf/*.cpp) $(wildcard GenotypeLikelihoods/*.cpp) $(wildcard PopulationTools/*.cpp) $(wildcard GLF/*.cpp) $(wildcard bamtools/api/*.cpp) $(wildcard bamtools/api/algorithms/*.cpp) $(wildcard bamtools/api/internal/bam/*.cpp) $(wildcard bamtools/api/internal/index/*.cpp) $(wildcard bamtools/api/internal/io/*.cpp) $(wildcard bamtools/api/internal/sam/*.cpp) $(wildcard bamtools/api/internal/utils/*.cpp) $(wildcard bamtools/utils/*.cpp)
GIT_HEADER = commonutilities/gitversion.cpp

OBJ = $(SRC:%.cpp=%.o)
BIN = atlas

.PHONY : all
all : $(BIN)

ifeq ($(ARM),)
BINFLAG = -lz -larmadillo
OBJFLAG = -std=c++1y
else
BINFLAG = -lz -lblas -llapack
OBJFLAG = -Iarmadillo-8.400.0/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -std=c++1y
endif

$(BIN): $(GIT_HEADER) $(OBJ)
	$(CXX) -O3 -o $(BIN) $(OBJ) $(BINFLAG)


$(GIT_HEADER): .git/HEAD .git/COMMIT_EDITMSG .git/index
	echo "#include \"gitversion.h\"" > $@
	echo "std::string getGitVersion(){" >> $@
	echo "return \"$(shell git rev-parse HEAD)\";" >> $@
	echo "}" >> $@

.git/COMMIT_EDITMSG :
	touch $@


%.o: %.cpp
	$(CXX) -O3 -c -I. -Ibamtools -Icommonutilities -IGLF -ITests -ISimulations -IVcf -IGenotypeLikelihoods -IPopulationTools $(OBJFLAG)  $< -o $@


.PHONY : clean
clean:
	rm -rf $(BIN) $(OBJ)
