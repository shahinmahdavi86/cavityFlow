# Compiler and flags
CXX = g++
CXX_FLAGS = -g -Wall
CXX_LIBS = -lm

# Directories
SRC_DIR = .
OBJ_DIR = obj

# Source and object files
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

# Output executable
BIN_FILE = cavity

# Default target
all: $(BIN_FILE)

# Link the object files to create the executable
$(BIN_FILE): $(OBJ_FILES)
	$(CXX) $(CXX_FLAGS) -o $@ $(OBJ_FILES) $(CXX_LIBS)

# Compile source files into object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXX_FLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -rf $(BIN_FILE) $(OBJ_DIR) results_vtk/*.vtk
