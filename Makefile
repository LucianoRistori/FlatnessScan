# ===============================================================
# Makefile for FlatnessScan
# Author: Luciano Ristori
# Updated: October 2025
#
# Description:
#   Builds the FlatnessScan program using the shared common/ module.
#   Dependencies: ROOT framework (macOS build via clang++).
# ===============================================================

# ---- Compiler and Flags ----
CXX      := clang++
CXXFLAGS := -O2 -Wall -Wextra -Wno-cpp -stdlib=libc++ -mmacosx-version-min=13.0 \
            -std=c++17 -Wno-c++17-extensions -pthread -m64
INCLUDES := -I../common -I. -I/Applications/root/include
LDFLAGS  := -L/Applications/root/lib -Wl,-rpath,/Applications/root/lib
LIBS     := -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps \
            -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore \
            -lThread -lMultiProc -lROOTDataFrame -lpthread -lm -ldl

# ---- Targets ----
TARGET    := flatnessscan
OBJS      := FlatnessScan.o ../common/Points.o

# ---- Default Rule ----
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(OBJS) -o $@ $(LDFLAGS) $(LIBS)

# ---- Pattern Rule for Compilation ----
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

../common/Points.o: ../common/Points.cpp ../common/Points.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ../common/Points.cpp -o ../common/Points.o

# ---- Cleaning ----
clean:
	@echo "Cleaning up..."
	rm -f *.o ../common/Points.o $(TARGET)
	find . -name "*.dSYM" -type d -exec rm -rf {} +

# ---- Convenience ----
rebuild: clean all

.PHONY: all clean rebuild
