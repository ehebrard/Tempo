cmake_minimum_required(VERSION 3.10)
project(tempo)
include(FetchContent)

# Also requires one of: libbfd (gnu binutils), libdwarf, libdw (elfutils)
FetchContent_Declare(backward
        GIT_REPOSITORY https://github.com/bombela/backward-cpp
        GIT_TAG v1.6)
FetchContent_MakeAvailable(backward)


# iterator library for zip and enumerate
string(COMPARE EQUAL ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)
string(COMPARE LESS ${CMAKE_CXX_COMPILER_VERSION} "16" LESS_THAN_16)
# Clangs ranges implementation is broken until clang-16 => only limited support of c++20 view features
if (${IS_CLANG} AND ${LESS_THAN_16})
    # This is a stable commit that does not make use of the c++20 view interface => no fancy pipe operator but at least
    # clang won't complain
    set(ITERATORS_COMMIT pre_cpp20_ranges)
else ()
    set(ITERATORS_COMMIT master)
endif ()

FetchContent_Declare(iterators
        GIT_REPOSITORY https://github.com/Timmifixedit/IteratorTools
        GIT_TAG ${ITERATORS_COMMIT}
)
FetchContent_MakeAvailable(iterators)
include_directories(${iterators_SOURCE_DIR})


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -mtune=native")
string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
if (BUILD_TYPE STREQUAL "debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
    message("Building for debug")
elseif (BUILD_TYPE STREQUAL "release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
    message("Building for release")
elseif (BUILD_TYPE STREQUAL "profile")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
    message("Building for profiling")
else ()
    message(SEND_ERROR "Unknown build type ${CMAKE_BUILD_TYPE}")
endif ()

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/cpp/*.cpp ${CMAKE_SOURCE_DIR}/src/cpp/util/*.cpp ${CMAKE_SOURCE_DIR}/src/cpp/constraints/*.cpp)

message("Source files:")
foreach(SRC ${SOURCES})
    message(\t${SRC})    
endforeach()

include_directories(${CMAKE_SOURCE_DIR}/src/header/)
include_directories(${CMAKE_SOURCE_DIR}/tools/)

# Examples
file(GLOB_RECURSE EXECUTABLES ${CMAKE_SOURCE_DIR}/src/examples/*.cpp)
message("Generating targets for the following cpp files:")
foreach(EXEC ${EXECUTABLES})
    get_filename_component(NAME ${EXEC} NAME_WLE)
    set(ALL_SOURCES ${SOURCES})
    message(\t${EXEC}\ ->\ target:\ ${NAME})
    if (BUILD_TYPE STREQUAL "debug")
        add_executable(${NAME} ${EXEC} ${ALL_SOURCES} ${BACKWARD_ENABLE})
        add_backward(${NAME})
    else()
        add_executable(${NAME} ${EXEC} ${ALL_SOURCES})
    endif()
endforeach()


## Creating a library, following tutorial here: https://stackoverflow.com/a/45843676
#message("Generating target for library")
#add_library(tempo STATIC ${SOURCES})
#include(GNUInstallDirs)
#install(TARGETS tempo
#    LIBRARY DESTINATION ${CMAKE_INSTALL_TEMPO}
#    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})