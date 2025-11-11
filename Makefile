# ================================================================
# FlatnessScan Makefile
# Author: Luciano Ristori
# Description:
#   Analyzes measured surface points, fits a plane, computes deviations,
#   and generates ROOT histograms and plots.
# ================================================================

CXX       = clang++
CXXFLAGS  = -O2 -Wall -Wextra -Wno-cpp -std=c++17 -stdlib=libc++ -pthread -m64 -mmacosx-version-min=13.0

# Automatically query ROOT for include and library paths
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

INCLUDES   = -I../common -I.

LDFLAGS    = -stdlib=libc++ -pthread -lm -ldl

SRCS       = FlatnessScan.cpp ../common/Points.cpp
OBJS       = $(SRCS:.cpp=.o)
TARGET     = flatnessScan

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(ROOTLIBS) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)
	find . -name "*.dSYM" -type d -exec rm -rf {} +
