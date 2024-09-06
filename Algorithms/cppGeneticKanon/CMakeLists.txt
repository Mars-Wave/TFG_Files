cmake_minimum_required(VERSION 3.28)
project(cppGeneticKanon)

set(CMAKE_CXX_STANDARD 20)

# Disable all optimizations
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

# Static linking for libstdc++
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")

add_executable(cppGeneticKanon
        geneticKAnonimity.cpp)