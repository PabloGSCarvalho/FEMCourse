# @file FemSC/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

project (FemSC)

add_subdirectory(headers)
add_subdirectory("src vazio")
add_subdirectory(tinyxml)

#Specify where the headers files must be found
include_directories(
	"${PROJECT_SOURCE_DIR}/headers"
	"${PROJECT_SOURCE_DIR}/tinyxml"
)

## inserting main definitions
#specifying the path to FemSC source code
set(FemSCSOURCE ${FemSC_SRCS} ${FemSC_tinyxml} ${FemSC_HEAD})

add_subdirectory(mainprograms)

## Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
endif()
