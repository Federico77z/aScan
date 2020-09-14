# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /usr/bin/cmake3

# The command to remove a file.
RM = /usr/bin/cmake3 -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/clap/aScan/git_repo/aScan/bamtools

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/clap/aScan/git_repo/aScan/bamtools/bin

# Include any dependencies generated for this target.
include src/toolkit/CMakeFiles/bamtools_cmd.dir/depend.make

# Include the progress variables for this target.
include src/toolkit/CMakeFiles/bamtools_cmd.dir/progress.make

# Include the compile flags for this target's objects.
include src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o: ../src/toolkit/bamtools_convert.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_convert.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_convert.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_convert.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o: ../src/toolkit/bamtools_count.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_count.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_count.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_count.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o: ../src/toolkit/bamtools_coverage.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_coverage.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_coverage.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_coverage.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o: ../src/toolkit/bamtools_filter.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_filter.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_filter.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_filter.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o: ../src/toolkit/bamtools_header.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_header.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_header.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_header.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o: ../src/toolkit/bamtools_index.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_index.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_index.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_index.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o: ../src/toolkit/bamtools_merge.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_merge.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_merge.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_merge.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o: ../src/toolkit/bamtools_random.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_random.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_random.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_random.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o: ../src/toolkit/bamtools_resolve.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_resolve.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_resolve.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_resolve.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o: ../src/toolkit/bamtools_revert.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_revert.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_revert.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_revert.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o: ../src/toolkit/bamtools_sort.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_sort.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_sort.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_sort.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o: ../src/toolkit/bamtools_split.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_split.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_split.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_split.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o: ../src/toolkit/bamtools_stats.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_stats.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_stats.cpp > CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools_stats.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.s

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o: src/toolkit/CMakeFiles/bamtools_cmd.dir/flags.make
src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o: ../src/toolkit/bamtools.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o -c /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools.cpp

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/bamtools_cmd.dir/bamtools.cpp.i"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools.cpp > CMakeFiles/bamtools_cmd.dir/bamtools.cpp.i

src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/bamtools_cmd.dir/bamtools.cpp.s"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit/bamtools.cpp -o CMakeFiles/bamtools_cmd.dir/bamtools.cpp.s

# Object files for target bamtools_cmd
bamtools_cmd_OBJECTS = \
"CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o" \
"CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o"

# External object files for target bamtools_cmd
bamtools_cmd_EXTERNAL_OBJECTS =

src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_convert.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_count.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_coverage.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_filter.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_header.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_index.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_merge.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_random.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_resolve.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_revert.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_sort.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_split.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools_stats.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/bamtools.cpp.o
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/build.make
src/toolkit/bamtools: src/api/libbamtools.a
src/toolkit/bamtools: src/utils/libbamtools-utils.a
src/toolkit/bamtools: src/third_party/jsoncpp/libjsoncpp.a
src/toolkit/bamtools: src/api/libbamtools.a
src/toolkit/bamtools: /usr/lib64/libz.so
src/toolkit/bamtools: src/toolkit/CMakeFiles/bamtools_cmd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/clap/aScan/git_repo/aScan/bamtools/bin/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Linking CXX executable bamtools"
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/bamtools_cmd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/toolkit/CMakeFiles/bamtools_cmd.dir/build: src/toolkit/bamtools

.PHONY : src/toolkit/CMakeFiles/bamtools_cmd.dir/build

src/toolkit/CMakeFiles/bamtools_cmd.dir/clean:
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit && $(CMAKE_COMMAND) -P CMakeFiles/bamtools_cmd.dir/cmake_clean.cmake
.PHONY : src/toolkit/CMakeFiles/bamtools_cmd.dir/clean

src/toolkit/CMakeFiles/bamtools_cmd.dir/depend:
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/clap/aScan/git_repo/aScan/bamtools /home/clap/aScan/git_repo/aScan/bamtools/src/toolkit /home/clap/aScan/git_repo/aScan/bamtools/bin /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit /home/clap/aScan/git_repo/aScan/bamtools/bin/src/toolkit/CMakeFiles/bamtools_cmd.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/toolkit/CMakeFiles/bamtools_cmd.dir/depend

