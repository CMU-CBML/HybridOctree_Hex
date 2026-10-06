cmake_minimum_required(VERSION 3.16)
project(HexGen VERSION 1.0)

# Mesh generation is compute intensive. Preserve an explicitly selected build
# type, and let multi-configuration generators select it at build time.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

set(SOURCES
    Main.cpp
    HexGen.cpp
    HexGen.h
    Initialization.h
    Mesh.cpp
    Mesh.h
    StaticVars.h
)

add_executable(HexGen ${SOURCES})
target_compile_features(HexGen PRIVATE cxx_std_11)
set_target_properties(HexGen PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)

# NATIVE preserves the original compiler's overload resolution. FLOATING matches
# a legacy build where unqualified abs(double) selected the floating overload.
set(HEXGEN_ABS_MODE "NATIVE" CACHE STRING "Legacy abs overload: NATIVE or FLOATING")
set_property(CACHE HEXGEN_ABS_MODE PROPERTY STRINGS NATIVE FLOATING)
if (HEXGEN_ABS_MODE STREQUAL "FLOATING")
    target_compile_definitions(HexGen PRIVATE HEXGEN_FLOATING_ABS)
elseif (NOT HEXGEN_ABS_MODE STREQUAL "NATIVE")
    message(FATAL_ERROR "HEXGEN_ABS_MODE must be NATIVE or FLOATING")
endif()

if (CMAKE_GENERATOR MATCHES "Visual Studio")
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT HexGen)
endif()

include(CTest)
if (BUILD_TESTING AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME legacy_equivalence COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../tools/verify_legacy.py"
        --suite fixtures --iterations 1000 --abs-mode "${HEXGEN_ABS_MODE}"
        --compiler "${CMAKE_CXX_COMPILER}"
        --output "${CMAKE_CURRENT_BINARY_DIR}/legacy-equivalence")
    set_tests_properties(legacy_equivalence PROPERTIES TIMEOUT 240)
elseif (BUILD_TESTING)
    message(STATUS "Legacy differential runner requires a Linux GCC/Clang environment")
endif()
