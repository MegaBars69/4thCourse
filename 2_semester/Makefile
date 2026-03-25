# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -g
INCLUDES = -I/usr/include/eigen3 -I.

# Source files
SRCS = main.cpp matrix.cpp mesh.cpp funcs.cpp P_gas.cpp P_scheme.cpp residuals.cpp table_tex_printer.cpp bicstab.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = solver

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Run the program with default parameters (adjust as needed)
run: $(TARGET)
	./$(TARGET) 0.1 1.0 1.0 1.0 1.0 50 50 10 1

# Debug build
debug: CXXFLAGS += -g -O0
debug: clean $(TARGET)

.PHONY: all clean run debug