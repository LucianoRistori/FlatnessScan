#----------------------------------------------
#  Makefile for FlatnessScan
#  Author: Luciano Ristori
#----------------------------------------------

# Compiler and flags
CXX      := clang++
CXXFLAGS := -std=c++17 -O2 -Wall -mmacosx-version-min=13.0
ROOTFLAGS := $(shell root-config --cflags --libs)

# Sources and executable
SOURCES  := flatnessScan.cpp Points.cpp
TARGET   := flatnessScan

#----------------------------------------------
all: $(TARGET)

$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(ROOTFLAGS)

clean:
	rm -f $(TARGET) *.o

#----------------------------------------------
# Usage:
#   make        -> build the program
#   make clean  -> remove executable and object files
#----------------------------------------------
