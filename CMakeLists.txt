cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(qn_discriminator)
set(CMAKE_CXX_FLAGS "-Wall")
# set (QnDiscriminator_VERSION_MAJOR 0)
# set (QnDiscriminator_VERSION_MINOR 9)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE RELEASE)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package( Boost REQUIRED COMPONENTS program_options )
find_package(ROOT REQUIRED)

if (NOT ROOT_FOUND)
    message(WARNING "ROOT was not found with CMAKE_PREFIX_PATH...")
    list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
    find_package(ROOT REQUIRED)
endif ()

if (ROOT_USE_FILE)
    message("-- ROOT_USE_FILE ${ROOT_USE_FILE}")
    include(${ROOT_USE_FILE})
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{QnTools})

find_package(QnTools REQUIRED)

#add_executable(draw src/main.cc src/draw_functions.cc  src/json_reader.cc src/file_manager.cc)
#target_link_libraries(draw QnTools::Base ${ROOT_LIBRARIES} ${Boost_LIBRARIES} )
#
set(SOURCES ShapeFitter.cpp
            Fitter.cpp
            GraphExtractor.cpp)
            
string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS "Helper.hpp")

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}
        ${Base_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${QnTools_INCLUDE_DIR}
)
        
set(PCM_FILE_NAME libQnDiscriminator)
message(STATUS "QnTools include ${QnTools_INCLUDE_DIR}")
ROOT_GENERATE_DICTIONARY( G__QnDiscriminator ${HEADERS}
        LINKDEF
        QnDiscriminatorLinkDef.h
        OPTIONS
#        -s ${PCM_FILE_NAME}
        -I ${QnTools_INCLUDE_DIR}/QnTools
        -I ${CMAKE_SOURCE_DIR}
        )
add_library(QnDiscriminator SHARED ${SOURCES} G__QnDiscriminator.cxx )
target_link_libraries(QnDiscriminator PRIVATE QnTools::Base ${ROOT_LIBRARIES})

add_executable(run run.cpp)
add_dependencies(run QnDiscriminator)
target_link_libraries(run QnTools::Base ${ROOT_LIBRARIES} QnDiscriminator)
