#==============================================================================
# Makefile for FlatnessScan
#==============================================================================

# Compiler
CXX = clang++

# Base compiler flags
CXXFLAGS = -O2 -Wall -Wextra -Wno-cpp \
           -std=c++17 -stdlib=libc++ \
           -mmacosx-version-min=13.0 \
           -pthread -m64

# ROOT configuration
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

# Common module
COMMON_DIR = ../common
INCLUDES   = -I$(COMMON_DIR)

# Sources
SRCS = flatnessScan.cpp \
       $(COMMON_DIR)/Points.cpp

# Objects
OBJS = $(SRCS:.cpp=.o)

# Target
TARGET = flatnessScan

# Where the built executable gets installed after every build, so it's
# always run from there (e.g. from a Dropbox working directory) rather
# than from this source directory.
BINDIR = $(HOME)/bin

# Default rule: build, then install into BINDIR. "all" is declared .PHONY
# below, so the install step runs on every "make" invocation, even when
# $(TARGET) itself was already up to date.
all: $(TARGET)
	@mkdir -p $(BINDIR)
	rm -f $(BINDIR)/$(TARGET)
	cp $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed $(TARGET) -> $(BINDIR)/$(TARGET)"

# Link
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(ROOTLIBS) -o $(TARGET)

# Compile
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
	find . -name "*.dSYM" -type d -exec rm -rf {} +

.PHONY: all clean
