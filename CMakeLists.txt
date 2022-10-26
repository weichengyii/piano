cmake_minimum_required(VERSION 3.24)
project(piano)

set(CMAKE_CXX_STANDARD 23)

add_executable(piano
        dwgs.cpp
        dwgs.h
        filter.cpp
        filter.h
        piano.cpp
        piano.h
        reverb.cpp
        reverb.h
        sndfile.cpp
        sndfile.h
        types.h
        hammer.cpp
        hammer.h
        )
