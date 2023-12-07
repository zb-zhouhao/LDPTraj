cmake_minimum_required(VERSION 3.17)
project(NNGram)
set(CMAKE_CXX_STANDARD_LIBRARIES "-lstdc++")
set(CMAKE_CXX_STANDARD 17)


# 如果想要生成的可执行文件拥有符号表，可以gdb调试，就直接加上这句
add_definitions("-Wall -g")
add_definitions("-lm")

# 加入ortools位置
set(CMAKE_PREFIX_PATH "./ortools")
find_package(ortools REQUIRED CONFIG)
include_directories(${PROJECT_SOURCE_DIR})

#set(BOOST_INCLUDEDIR "/usr/local/include/boost")
#set(BOOST_LIBDIR "/usr/local/lib")

find_package(Boost 1.77 REQUIRED COMPONENTS math_c99)
#find_package(Boost REQUIRED COMPONENTS math_c99)
#find_package(Boost)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})

    MESSAGE( STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    MESSAGE( STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    MESSAGE( STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
endif()

add_executable(NNGram main.cpp util.h io/FileReader.cpp io/FileReader.h Parameters.cpp Parameters.h io/LogFile.cpp io/LogFile.h  NG/NGramSolver.cpp NG/NGramSolver.h
        Triple.h GI/GeoISolver.cpp GI/GeoISolver.h io/FileWriter.cpp io/FileWriter.h Local/ExpMechSolver.cpp Local/ExpMechSolver.h KP/KeyPointSolver.cpp KP/KeyPointSolver.h)

target_link_libraries(${PROJECT_NAME} PRIVATE ortools::ortools Boost::math_c99 -lm)
#target_link_libraries(${PROJECT_NAME} PRIVATE  Boost::math_c99)
#target_link_libraries(${PROJECT_NAME} )
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} -lm)
