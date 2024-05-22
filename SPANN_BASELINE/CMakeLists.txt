cmake_minimum_required(VERSION 3.21)
project(SPANN_BASELINE)

set(CMAKE_CXX_STANDARD 14)

add_executable(SPANN_BASELINE main.cpp K-means.h Point.h Baseline.h monitor.h)
