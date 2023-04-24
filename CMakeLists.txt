cmake_minimum_required (VERSION 2.6)
project (APPSPAM)

#SET(GCC_COMPILE_FLAGS "-O3 -std=c++11 -g -Wno-write-strings -fopenmp")
#add_definitions(${GCC_COMPILE_FLAGS})
set(CMAKE_CXX_STANDARD 17)  # enable C++17 standard

# include
include_directories(./include ./src)

# target
add_executable(appspam
        src/Algorithms.cpp include/Algorithms.h
        src/APPSPAM.cpp
        src/BucketManager.cpp include/BucketManager.h
        src/Crc32.cpp include/Crc32.h
        src/GenomeManager.cpp include/GenomeManager.h
        src/GlobalParameters.cpp include/GlobalParameters.h
        src/Match.cpp include/Match.h
        src/MatchManager.cpp include/MatchManager.h
        src/Node.cpp include/Node.h
        src/Pattern.cpp include/Pattern.h
        src/Placement.cpp include/Placement.h
        src/ReadManager.cpp include/ReadManager.h
        src/Scoring.cpp include/Scoring.h
        src/Seed.cpp include/Seed.h
        src/SeqIO.cpp include/SeqIO.h
        src/Sequence.cpp include/Sequence.h
        src/SubstitutionMatrix.h
        src/Tree.cpp include/Tree.h)

# Enable all warnings, let's look at them <3
target_compile_options(appspam PRIVATE -Wall -Wextra)


# OpenMP
FIND_PACKAGE(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
include_directories(SYSTEM ${OpenMP_INCLUDE_PATH})
