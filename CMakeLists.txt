cmake_minimum_required(VERSION 3.25)

project(IMASH LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)

set(ENV{PKG_CONFIG_PATH}
    "$ENV{HOME}/IMAS-Cpp/installdir/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}"
)

pkg_check_modules(IMASCPP REQUIRED IMPORTED_TARGET al-cpp)

add_library(IMASH STATIC
    src/IMASH.cpp
    src/TokamakData.cpp
)

target_include_directories(IMASH PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(IMASH PUBLIC
    PkgConfig::IMASCPP
    Threads::Threads
)

add_executable(imash_test
    src/main.cpp
)
target_compile_options(IMASH PUBLIC -pthread)
target_link_options(IMASH PUBLIC -pthread)
target_link_libraries(imash_test PRIVATE IMASH)
set_target_properties(imash_test PROPERTIES
    BUILD_RPATH "${IMASCPP_LIBRARY_DIRS}"
    INSTALL_RPATH "${IMASCPP_LIBRARY_DIRS}"
    SKIP_BUILD_RPATH FALSE
)

message(STATUS "IMASCPP_LIBRARY_DIRS = ${IMASCPP_LIBRARY_DIRS}")
message(STATUS "IMASCPP_LIBRARIES    = ${IMASCPP_LIBRARIES}")