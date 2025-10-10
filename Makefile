# =========================
# Makefile for makeHists
# =========================

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -mmacosx-version-min=13.0

# ROOT flags
ROOTCFLAGS = `root-config --cflags`
ROOTLIBS   = `root-config --libs`

# Source files
SRCS = flatnessScan.cpp Points.cpp
OBJS = $(SRCS:.cpp=.o)

# Target executable
TARGET = flatnessScan

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(ROOTLIBS)

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< $(ROOTCFLAGS)

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
