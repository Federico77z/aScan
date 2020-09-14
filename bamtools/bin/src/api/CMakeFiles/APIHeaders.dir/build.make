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

# Utility rule file for APIHeaders.

# Include the progress variables for this target.
include src/api/CMakeFiles/APIHeaders.dir/progress.make

APIHeaders: src/api/CMakeFiles/APIHeaders.dir/build.make
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/api_global.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/api_global.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamAlgorithms.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamAlgorithms.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamAlignment.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamAlignment.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamAux.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamAux.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamConstants.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamConstants.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamIndex.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamIndex.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamMultiReader.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamMultiReader.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamReader.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamReader.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/BamWriter.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/BamWriter.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/IBamIODevice.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/IBamIODevice.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamConstants.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamConstants.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamHeader.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamHeader.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamProgram.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamProgram.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamProgramChain.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamProgramChain.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamReadGroup.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamReadGroup.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamReadGroupDictionary.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamReadGroupDictionary.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamSequence.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamSequence.h
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && /usr/bin/cmake3 -E copy_if_different /home/clap/aScan/git_repo/aScan/bamtools/src/api/SamSequenceDictionary.h /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/include/api/SamSequenceDictionary.h
.PHONY : APIHeaders

# Rule to build all files generated by this target.
src/api/CMakeFiles/APIHeaders.dir/build: APIHeaders

.PHONY : src/api/CMakeFiles/APIHeaders.dir/build

src/api/CMakeFiles/APIHeaders.dir/clean:
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api && $(CMAKE_COMMAND) -P CMakeFiles/APIHeaders.dir/cmake_clean.cmake
.PHONY : src/api/CMakeFiles/APIHeaders.dir/clean

src/api/CMakeFiles/APIHeaders.dir/depend:
	cd /home/clap/aScan/git_repo/aScan/bamtools/bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/clap/aScan/git_repo/aScan/bamtools /home/clap/aScan/git_repo/aScan/bamtools/src/api /home/clap/aScan/git_repo/aScan/bamtools/bin /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api /home/clap/aScan/git_repo/aScan/bamtools/bin/src/api/CMakeFiles/APIHeaders.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/api/CMakeFiles/APIHeaders.dir/depend
