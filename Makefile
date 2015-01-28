# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:
.PHONY : .NOTPARALLEL

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /auto-home/stud/me31kove/Documents/siwir_git/siwir_ex04

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /auto-home/stud/me31kove/Documents/siwir_git/siwir_ex04

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /auto-home/stud/me31kove/Documents/siwir_git/siwir_ex04/CMakeFiles /auto-home/stud/me31kove/Documents/siwir_git/siwir_ex04/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /auto-home/stud/me31kove/Documents/siwir_git/siwir_ex04/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named heat

# Build rule for target.
heat: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 heat
.PHONY : heat

# fast build rule for target.
heat/fast:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/build
.PHONY : heat/fast

src/Array.o: src/Array.cc.o
.PHONY : src/Array.o

# target to build an object file
src/Array.cc.o:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Array.cc.o
.PHONY : src/Array.cc.o

src/Array.i: src/Array.cc.i
.PHONY : src/Array.i

# target to preprocess a source file
src/Array.cc.i:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Array.cc.i
.PHONY : src/Array.cc.i

src/Array.s: src/Array.cc.s
.PHONY : src/Array.s

# target to generate assembly for a file
src/Array.cc.s:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Array.cc.s
.PHONY : src/Array.cc.s

src/CGSolver.o: src/CGSolver.cc.o
.PHONY : src/CGSolver.o

# target to build an object file
src/CGSolver.cc.o:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/CGSolver.cc.o
.PHONY : src/CGSolver.cc.o

src/CGSolver.i: src/CGSolver.cc.i
.PHONY : src/CGSolver.i

# target to preprocess a source file
src/CGSolver.cc.i:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/CGSolver.cc.i
.PHONY : src/CGSolver.cc.i

src/CGSolver.s: src/CGSolver.cc.s
.PHONY : src/CGSolver.s

# target to generate assembly for a file
src/CGSolver.cc.s:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/CGSolver.cc.s
.PHONY : src/CGSolver.cc.s

src/Debug.o: src/Debug.cc.o
.PHONY : src/Debug.o

# target to build an object file
src/Debug.cc.o:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Debug.cc.o
.PHONY : src/Debug.cc.o

src/Debug.i: src/Debug.cc.i
.PHONY : src/Debug.i

# target to preprocess a source file
src/Debug.cc.i:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Debug.cc.i
.PHONY : src/Debug.cc.i

src/Debug.s: src/Debug.cc.s
.PHONY : src/Debug.s

# target to generate assembly for a file
src/Debug.cc.s:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/Debug.cc.s
.PHONY : src/Debug.cc.s

src/HeatSolver.o: src/HeatSolver.cc.o
.PHONY : src/HeatSolver.o

# target to build an object file
src/HeatSolver.cc.o:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/HeatSolver.cc.o
.PHONY : src/HeatSolver.cc.o

src/HeatSolver.i: src/HeatSolver.cc.i
.PHONY : src/HeatSolver.i

# target to preprocess a source file
src/HeatSolver.cc.i:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/HeatSolver.cc.i
.PHONY : src/HeatSolver.cc.i

src/HeatSolver.s: src/HeatSolver.cc.s
.PHONY : src/HeatSolver.s

# target to generate assembly for a file
src/HeatSolver.cc.s:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/HeatSolver.cc.s
.PHONY : src/HeatSolver.cc.s

src/main.o: src/main.cc.o
.PHONY : src/main.o

# target to build an object file
src/main.cc.o:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/main.cc.o
.PHONY : src/main.cc.o

src/main.i: src/main.cc.i
.PHONY : src/main.i

# target to preprocess a source file
src/main.cc.i:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/main.cc.i
.PHONY : src/main.cc.i

src/main.s: src/main.cc.s
.PHONY : src/main.s

# target to generate assembly for a file
src/main.cc.s:
	$(MAKE) -f CMakeFiles/heat.dir/build.make CMakeFiles/heat.dir/src/main.cc.s
.PHONY : src/main.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... heat"
	@echo "... rebuild_cache"
	@echo "... src/Array.o"
	@echo "... src/Array.i"
	@echo "... src/Array.s"
	@echo "... src/CGSolver.o"
	@echo "... src/CGSolver.i"
	@echo "... src/CGSolver.s"
	@echo "... src/Debug.o"
	@echo "... src/Debug.i"
	@echo "... src/Debug.s"
	@echo "... src/HeatSolver.o"
	@echo "... src/HeatSolver.i"
	@echo "... src/HeatSolver.s"
	@echo "... src/main.o"
	@echo "... src/main.i"
	@echo "... src/main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

