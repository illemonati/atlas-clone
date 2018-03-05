#make file for atlas

SRC = $(wildcard *.cpp) $(wildcard *.C) $(wildcard bamtools/api/*.cpp) $(wildcard bamtools/api/algorithms/*.cpp) $(wildcard bamtools/api/internal/bam/*.cpp) $(wildcard bamtools/api/internal/index/*.cpp) $(wildcard bamtools/api/internal/io/*.cpp) $(wildcard bamtools/api/internal/sam/*.cpp) $(wildcard bamtools/api/internal/utils/*.cpp) $(wildcard bamtools/utils/*.cpp) 

OBJ = $(SRC:%.cpp=%.o)

BIN = atlas

all:	$(BIN)

#replace command below by the following if armadillo cannot be installed system-wide or you get linking errors:
#$(CXX) -O3 -o $(BIN) $(OBJ) -lz -lblas -llapack

$(BIN):	$(OBJ)
	$(CXX) -O3 -o $(BIN) $(OBJ) -lz -larmadillo
	

#replace command below by the following if armadillo cannot be installed system-wide or you get linking errors:
#$(CXX) -O3 -c -Ibamtools -Iarmadillo-VERSION/include -DARMA_DONT_USE_WRAPPER -lblas -llapack -std=c++1y $< -o $@

%.o: %.cpp
	$(CXX) -O3 -c -Ibamtools -std=c++1y $< -o $@
	

clean:
	rm -rf $(BIN) $(OBJ)