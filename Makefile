#============================================================
# Makefile for FlatnessScan and GridFinder tools
#============================================================

#------------------------------------------------------------
# Compiler and ROOT configuration
#------------------------------------------------------------
CXX      := clang++
CXXFLAGS := -O2 -Wall -Wextra -Wno-cpp -stdlib=libc++ -mmacosx-version-min=13.0 -std=c++17 
CXXFLAGS += -Wno-c++17-extensions -I../common
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
LDFLAGS := -pthread -lm -ldl -Wl,-rpath,/Applications/root/lib -std=c++17

# Include directories
INCLUDES := -I. -I/Applications/root/include

# Source files
SRC := FlatnessScan.cpp ../common/Points.cpp
OBJ := $(SRC:.cpp=.o)

# Executables
EXE := flatnessscan
TEST_EXE := testGridFinder

#============================================================
# Default target: build both main and test executables
#============================================================
all: $(EXE) $(TEST_EXE)

#------------------------------------------------------------
# Main program: flatnessScan
#------------------------------------------------------------
$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $@ $(ROOTLIBS) $(LDFLAGS)

#------------------------------------------------------------
# Test program for GridFinder
#------------------------------------------------------------
$(TEST_EXE): testGridFinder.cpp Points.o GridFinder.h
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) -o $@ testGridFinder.cpp Points.o $(ROOTLIBS) $(LDFLAGS)

#------------------------------------------------------------
# Object file build rule
#------------------------------------------------------------
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) -c $< -o $@

#------------------------------------------------------------
# Cleaning
#------------------------------------------------------------
clean:
	rm -f $(OBJ) $(EXE) $(TEST_EXE)
	find . -name "*.dSYM" -type d -exec rm -rf {} +

#------------------------------------------------------------
# Notes
#------------------------------------------------------------
# - GridFinder.h and version.h are header-only; no .o file required
# - ROOT must be initialized in your shell (source initRoot)
# - Run:
#       make            # builds both executables
#       make flatnessScan
#       make testGridFinder
#       make clean
#============================================================

