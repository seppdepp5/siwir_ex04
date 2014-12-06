CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

ANSI_CFLAGS  = -std=c++0x

CFLAGS   = -O3  $(ANSI_CFLAGS) -DNDEBUG -Wall -Wextra -Wshadow -Werror -fopenmp -march=native
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
CPPFLAGS = -std=c++0x
LFLAGS   = -fopenmp
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
