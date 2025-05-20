# Makefile for building the Robust Savitzky-Golay Filter library
# Author: Mapoet
# Date: 2025-05-20

# Compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2018 -fimplicit-none -fcheck=all -fbacktrace
FFLAGS_RELEASE = -O3 -Wall -std=f2018 -fimplicit-none
INCLUDES = -I$(INCLUDE_DIR)

# Library dependencies
LIBS = -lopenblas -llapack

# Directories
INCLUDE_DIR = include
SRC_DIR = src
BUILD_DIR = build
LIB_DIR = lib

# Source files
SRC_FILES = $(SRC_DIR)/sg_coeff.f90 \
            $(SRC_DIR)/sg_filter.f90

INCLUDE_FILES = $(INCLUDE_DIR)/sglib_mod.f90

# Object files
OBJ_FILES = $(BUILD_DIR)/sglib_mod.o \
            $(BUILD_DIR)/sg_coeff.o \
            $(BUILD_DIR)/sg_filter.o

# Library output
LIB_STATIC = $(LIB_DIR)/libsglib.a
LIB_SHARED = $(LIB_DIR)/libsglib.so

# Test program
TEST_PROG = $(BUILD_DIR)/test_sg_filter
TEST_CUSTOM = $(BUILD_DIR)/test_sg_filter_data
# Targets
.PHONY: all clean test debug release

all: debug

debug: FFLAGS += -g
debug: directories $(LIB_STATIC) $(LIB_SHARED) $(TEST_PROG)

release: FFLAGS = $(FFLAGS_RELEASE)
release: directories $(LIB_STATIC) $(LIB_SHARED) $(TEST_PROG)

# Create required directories
directories:
	mkdir -p $(BUILD_DIR) $(LIB_DIR)

# Static library build
$(LIB_STATIC): $(OBJ_FILES)
	ar rcs $@ $^

# Shared library build
$(LIB_SHARED): $(OBJ_FILES)
	$(FC) -shared -o $@ $^ $(LIBS)

# Test program build
$(TEST_PROG): $(BUILD_DIR)/main_test.o $(LIB_STATIC)
	$(FC) -o $@ $< $(LIB_STATIC) $(LIBS)
$(TEST_CUSTOM): $(BUILD_DIR)/test_sg_filter_data.o $(LIB_STATIC)
	$(FC) -o $@ $< $(LIB_STATIC) $(LIBS)

# Compile module file
$(BUILD_DIR)/sglib_mod.o: $(INCLUDE_FILES)
	$(FC) $(FFLAGS) -c $< -o $@

# Compile submodule files
$(BUILD_DIR)/sg_coeff.o: $(SRC_DIR)/sg_coeff.f90 $(BUILD_DIR)/sglib_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDES)

$(BUILD_DIR)/sg_filter.o: $(SRC_DIR)/sg_filter.f90 $(BUILD_DIR)/sglib_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDES)

# Compile test program
$(BUILD_DIR)/main_test.o: $(SRC_DIR)/main_test.f90 $(BUILD_DIR)/sglib_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDES)

$(BUILD_DIR)/test_sg_filter_data.o: $(SRC_DIR)/test_sg_filter_data.f90 $(BUILD_DIR)/sglib_mod.o
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDES)

# Run test
test: $(TEST_PROG)
	cd $(BUILD_DIR) && ./test_sg_filter

# Run test with test data
test_data: $(TEST_CUSTOM)
	cd $(BUILD_DIR) && ./test_sg_filter_data -i ./test_data.dat -r 2.0 -p 3 -w 150

# Clean up
clean:
	rm -rf $(BUILD_DIR) $(LIB_DIR) sg_filter_results.dat 