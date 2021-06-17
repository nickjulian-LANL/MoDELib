# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build

# Include any dependencies generated for this target.
include Lattices/CMakeFiles/MoDELibLattices.dir/depend.make

# Include the progress variables for this target.
include Lattices/CMakeFiles/MoDELibLattices.dir/progress.make

# Include the compile flags for this target's objects.
include Lattices/CMakeFiles/MoDELibLattices.dir/flags.make

Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o: Lattices/CMakeFiles/MoDELibLattices.dir/flags.make
Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o: /Users/giacomo/Documents/MODELIB2/src/Lattices/Lattice.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o"
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && /opt/local/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o -c /Users/giacomo/Documents/MODELIB2/src/Lattices/Lattice.cpp

Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MoDELibLattices.dir/Lattice.cpp.i"
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && /opt/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giacomo/Documents/MODELIB2/src/Lattices/Lattice.cpp > CMakeFiles/MoDELibLattices.dir/Lattice.cpp.i

Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MoDELibLattices.dir/Lattice.cpp.s"
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && /opt/local/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giacomo/Documents/MODELIB2/src/Lattices/Lattice.cpp -o CMakeFiles/MoDELibLattices.dir/Lattice.cpp.s

# Object files for target MoDELibLattices
MoDELibLattices_OBJECTS = \
"CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o"

# External object files for target MoDELibLattices
MoDELibLattices_EXTERNAL_OBJECTS =

Lattices/libMoDELibLattices.a: Lattices/CMakeFiles/MoDELibLattices.dir/Lattice.cpp.o
Lattices/libMoDELibLattices.a: Lattices/CMakeFiles/MoDELibLattices.dir/build.make
Lattices/libMoDELibLattices.a: Lattices/CMakeFiles/MoDELibLattices.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libMoDELibLattices.a"
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && $(CMAKE_COMMAND) -P CMakeFiles/MoDELibLattices.dir/cmake_clean_target.cmake
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MoDELibLattices.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Lattices/CMakeFiles/MoDELibLattices.dir/build: Lattices/libMoDELibLattices.a

.PHONY : Lattices/CMakeFiles/MoDELibLattices.dir/build

Lattices/CMakeFiles/MoDELibLattices.dir/clean:
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices && $(CMAKE_COMMAND) -P CMakeFiles/MoDELibLattices.dir/cmake_clean.cmake
.PHONY : Lattices/CMakeFiles/MoDELibLattices.dir/clean

Lattices/CMakeFiles/MoDELibLattices.dir/depend:
	cd /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator /Users/giacomo/Documents/MODELIB2/src/Lattices /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices /Users/giacomo/Documents/MODELIB2/tools/MicrostructureGenerator/build/Lattices/CMakeFiles/MoDELibLattices.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Lattices/CMakeFiles/MoDELibLattices.dir/depend

