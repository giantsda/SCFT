# Set the name of the project and target:
SET(TARGET "scft")
PROJECT(${TARGET}) #set project

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#  FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#  FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
#  SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
SET(TARGET_SRC
  scft.cc
  )

SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
SET_PROPERTY(GLOBAL PROPERTY CMAKE_USE_RELATIVE_PATHS TRUE) 

#file(GLOB SOURCES "src/*.cc")
set(SOURCES ../gaussj.c ../broydn.c ../nrutil.c) # I need them

set_source_files_properties( ${SOURCES} PROPERTIES LANGUAGE CXX ) #use g++ to bulid them

ADD_LIBRARY(mylib  ${SOURCES}) # for link to scft

 
# No sensible run target can be defined for this example step
SET(TARGET_RUN "./scft")

include_directories( ${CMAKE_CURRENT_BINARY_DIR}/../ ) #provide .h
 

find_library(.. libNR_CPP.so)
link_directories(${CMAKE_CURRENT_BINARY_DIR}/../) # -L..
 
 

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

FIND_PACKAGE(deal.II 8.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()


DEAL_II_INITIALIZE_CACHED_VARIABLES()
SET(CLEAN_UP_FILES *.o scft *.vtk solution_store.txt solution_yita_1D*.txt *.msh)
 
DEAL_II_INVOKE_AUTOPILOT()

TARGET_LINK_LIBRARIES(scft  libNR_CPP.so mylib)
 


