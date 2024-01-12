cmake_minimum_required(VERSION 3.27.4)
project(HexGen VERSION 1.0)

# Add all source files to the project
set(SOURCES
    Main.cpp
    CellQueue.h
    HexGen.cpp
    HexGen.h
    Initialization.h
    Mesh.cpp
    Mesh.h
    StaticVars.h
)

# Add executable and link source files
add_executable(HexGen ${SOURCES})
