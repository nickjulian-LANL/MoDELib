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
CMAKE_SOURCE_DIR = /Users/giacomo/Downloads/TRI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/giacomo/Downloads/TRI

# Include any dependencies generated for this target.
include CMakeFiles/tri.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tri.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tri.dir/flags.make

CMakeFiles/tri.dir/src/main.cpp.o: CMakeFiles/tri.dir/flags.make
CMakeFiles/tri.dir/src/main.cpp.o: src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/giacomo/Downloads/TRI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tri.dir/src/main.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tri.dir/src/main.cpp.o -c /Users/giacomo/Downloads/TRI/src/main.cpp

CMakeFiles/tri.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tri.dir/src/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giacomo/Downloads/TRI/src/main.cpp > CMakeFiles/tri.dir/src/main.cpp.i

CMakeFiles/tri.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tri.dir/src/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giacomo/Downloads/TRI/src/main.cpp -o CMakeFiles/tri.dir/src/main.cpp.s

CMakeFiles/tri.dir/src/assert.cpp.o: CMakeFiles/tri.dir/flags.make
CMakeFiles/tri.dir/src/assert.cpp.o: src/assert.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/giacomo/Downloads/TRI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/tri.dir/src/assert.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tri.dir/src/assert.cpp.o -c /Users/giacomo/Downloads/TRI/src/assert.cpp

CMakeFiles/tri.dir/src/assert.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tri.dir/src/assert.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giacomo/Downloads/TRI/src/assert.cpp > CMakeFiles/tri.dir/src/assert.cpp.i

CMakeFiles/tri.dir/src/assert.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tri.dir/src/assert.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giacomo/Downloads/TRI/src/assert.cpp -o CMakeFiles/tri.dir/src/assert.cpp.s

CMakeFiles/tri.dir/src/del_impl.cpp.o: CMakeFiles/tri.dir/flags.make
CMakeFiles/tri.dir/src/del_impl.cpp.o: src/del_impl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/giacomo/Downloads/TRI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/tri.dir/src/del_impl.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tri.dir/src/del_impl.cpp.o -c /Users/giacomo/Downloads/TRI/src/del_impl.cpp

CMakeFiles/tri.dir/src/del_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tri.dir/src/del_impl.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/giacomo/Downloads/TRI/src/del_impl.cpp > CMakeFiles/tri.dir/src/del_impl.cpp.i

CMakeFiles/tri.dir/src/del_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tri.dir/src/del_impl.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/giacomo/Downloads/TRI/src/del_impl.cpp -o CMakeFiles/tri.dir/src/del_impl.cpp.s

# Object files for target tri
tri_OBJECTS = \
"CMakeFiles/tri.dir/src/main.cpp.o" \
"CMakeFiles/tri.dir/src/assert.cpp.o" \
"CMakeFiles/tri.dir/src/del_impl.cpp.o"

# External object files for target tri
tri_EXTERNAL_OBJECTS =

tri: CMakeFiles/tri.dir/src/main.cpp.o
tri: CMakeFiles/tri.dir/src/assert.cpp.o
tri: CMakeFiles/tri.dir/src/del_impl.cpp.o
tri: CMakeFiles/tri.dir/build.make
tri: CMakeFiles/tri.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/giacomo/Downloads/TRI/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable tri"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tri.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tri.dir/build: tri

.PHONY : CMakeFiles/tri.dir/build

CMakeFiles/tri.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tri.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tri.dir/clean

CMakeFiles/tri.dir/depend:
	cd /Users/giacomo/Downloads/TRI && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/giacomo/Downloads/TRI /Users/giacomo/Downloads/TRI /Users/giacomo/Downloads/TRI /Users/giacomo/Downloads/TRI /Users/giacomo/Downloads/TRI/CMakeFiles/tri.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tri.dir/depend

