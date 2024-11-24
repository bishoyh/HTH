
# Makefile for compiling main.cc

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3 -static -std=c++11 -pthread

# Target executable
TARGET = hth

# Source files
SOURCES = main.cc

# Object files
OBJECTS = $(SOURCES:.cc=.o)

# Build the target executable
$(TARGET): $(OBJECTS)
    $(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Clean up generated files
clean:
    rm -f $(TARGET) $(OBJECTS)

# Compile source files into object files
%.o: %.cc
    $(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
