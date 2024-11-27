cmake_minimum_required(VERSION 3.10)
project(tempo)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
include(GetGitRevisionDescription)
include(FetchContent)

# Also requires one of: libbfd (gnu binutils), libdwarf, libdw (elfutils)
FetchContent_Declare(backward
        GIT_REPOSITORY https://github.com/bombela/backward-cpp
        GIT_TAG v1.6)
FetchContent_MakeAvailable(backward)

# Json lib
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)

# iterator library for zip and enumerate
string(COMPARE EQUAL ${CMAKE_CXX_COMPILER_ID} "Clang" IS_CLANG)
string(COMPARE LESS ${CMAKE_CXX_COMPILER_VERSION} "16" LESS_THAN_16)
# Clangs ranges implementation is broken until clang-16 => only limited support of c++20 view features
if (${IS_CLANG} AND ${LESS_THAN_16})
    add_compile_definitions(__STD_RANGES_DISABLED__)
endif ()

FetchContent_Declare(iterators
        GIT_REPOSITORY https://github.com/Timmifixedit/IteratorTools
        GIT_TAG master
)
FetchContent_MakeAvailable(iterators)
include_directories(${iterators_SOURCE_DIR})


find_package(TorchScatter QUIET)
find_package(TorchSparse QUIET)

option(ENABLE_TORCH "Build with torch. If set to off (default), you are unable to use features that depend on torch")
option(PRETTY_JSON "Save json files in a more human readable format at the cost of higher memory usage (default on)" ON)
option(CLUSTER "Build for cluster without architecture specific optimizations" OFF)
set(TORCH_ENABLED FALSE)
#set(PRETTY_JSON FALSE)

if (${ENABLE_TORCH})
    if (NOT TorchScatter_FOUND OR NOT TorchSparse_FOUND)
        message(SEND_ERROR "You requested the use of torch but the library was not found on the system")
    elseif (${ENABLE_TORCH})
        set(TORCH_ENABLED TRUE)
    endif ()
endif ()

if (${PRETTY_JSON})
    add_compile_definitions(__JSON_INDENT__=2)
else ()
    add_compile_definitions(__JSON_INDENT__=)
endif ()

add_compile_definitions(__TEST_DATA_DIR__="${CMAKE_SOURCE_DIR}/Tests/data")
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${CMAKE_SOURCE_DIR}/git_sha.cpp.in"
        "${CMAKE_CURRENT_BINARY_DIR}/git_sha.cpp" @ONLY)

if (${TORCH_ENABLED})
    add_compile_definitions(__TORCH_ENABLED__=1)
    set(TORCH_LIB TorchScatter::TorchScatter TorchSparse::TorchSparse)
    set(CMAKE_CXX_FLAGS ${TORCH_CXX_FLAGS})
else ()
    add_compile_definitions(__TORCH_ENABLED__=0)
endif ()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
if (${CLUSTER})
    message("Building for cluster")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
else ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mtune=native -march=native")
    if (BUILD_TYPE STREQUAL "debug")
        add_compile_definitions(__DEBUG_BUILD__)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -g")
        message("Building for debug")
    elseif (BUILD_TYPE STREQUAL "release")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
        message("Building for release")
    else ()
        message(SEND_ERROR "Unknown build type ${CMAKE_BUILD_TYPE}")
    endif ()
endif ()

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/cpp/*.cpp ${CMAKE_SOURCE_DIR}/src/cpp/util/*.cpp
        ${CMAKE_SOURCE_DIR}/src/cpp/constraints/*.cpp ${CMAKE_SOURCE_DIR}/src/cpp/heuristics/*.cpp
        ${CMAKE_SOURCE_DIR}/src/cpp/heuristics/LNS/*.cpp ${CMAKE_CURRENT_BINARY_DIR}/git_sha.cpp)

file(GLOB EXAMPLE_HELPER_SOURCES ${CMAKE_SOURCE_DIR}/src/examples/helpers/*.cpp)

file(GLOB TORCH_SOURCES ${CMAKE_SOURCE_DIR}/src/cpp/nn/*.cpp)

message("compiler flags: ${CMAKE_CXX_FLAGS}")

message("Source files:")
foreach(SRC ${SOURCES})
    message(\t${SRC})
endforeach()

include_directories(${CMAKE_SOURCE_DIR}/src/header/)
include_directories(${CMAKE_SOURCE_DIR}/tools/)

# Examples
file(GLOB EXECUTABLES ${CMAKE_SOURCE_DIR}/src/examples/*.cpp ${CMAKE_SOURCE_DIR}/src/examples/torch/*.cpp)
message("Generating targets for the following cpp files:")
foreach(EXEC ${EXECUTABLES})
    get_filename_component(NAME ${EXEC} NAME_WLE)
    string(FIND ${EXEC} "torch" TORCH_EXEC)
    set(LIB_SOURCES ${SOURCES})
    if (NOT ${TORCH_EXEC} EQUAL -1)
        if (NOT ${TORCH_ENABLED})
            message("--- Info: Torch features disabled. Skipping executable " ${EXEC} " ---")
            continue()
        endif ()
        set(LIB_SOURCES ${LIB_SOURCES} ${TORCH_SOURCES})
    endif ()
    message(\t${EXEC}\ ->\ target:\ ${NAME})
    if (BUILD_TYPE STREQUAL "debug")
        add_executable(${NAME} ${EXEC} ${LIB_SOURCES} ${EXAMPLE_HELPER_SOURCES} ${BACKWARD_ENABLE})
        add_backward(${NAME})
    else()
        add_executable(${NAME} ${EXEC} ${LIB_SOURCES} ${EXAMPLE_HELPER_SOURCES})
    endif()

    if (NOT ${TORCH_EXEC} EQUAL -1 AND ${TORCH_ENABLED})
        target_include_directories(${NAME} PRIVATE ${TorchScatter_INCLUDE_DIR})
        target_include_directories(${NAME} PRIVATE ${TorchSparse_INCLUDE_DIR})
        target_link_libraries(${NAME} PRIVATE ${TORCH_LIB})
    endif ()
    target_link_libraries(${NAME} PRIVATE nlohmann_json::nlohmann_json)
endforeach()

add_library(tempo STATIC ${SOURCES})
target_link_libraries(tempo PRIVATE nlohmann_json::nlohmann_json)


add_subdirectory(Tests)