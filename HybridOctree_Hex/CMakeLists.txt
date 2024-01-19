cmake_minimum_required(VERSION 3.27.4)
project(HexGen VERSION 1.0)

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

if (CMAKE_GENERATOR MATCHES "Visual Studio")
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT HexGen)
endif()