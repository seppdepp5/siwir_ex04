CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( heat )

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})


SET( SRCS

   	src/main.cc
 	src/HeatSolver.cc 
	src/CGSolver.cc 
	src/Array.cc 
	src/Debug.cc 
   # more *.cpp files follow here
)


SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -std=c++0x")

OPTION( ARCHITECTURE_OPTIMIZATION "enables/disables flag -march=native (ON by default)" ON )

OPTION( ENABLE_WERROR "enables/disables flag -Werror (ON by default)" ON )

IF( CMAKE_COMPILER_IS_GNUCC )
   IF( ARCHITECTURE_OPTIMIZATION )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native" )
   ENDIF()
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCC )
   IF( ENABLE_WERROR )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror" )
   ENDIF()
ENDIF()

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")

ADD_EXECUTABLE( heat ${SRCS} )

target_link_libraries(heat ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(heat PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(heat PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
