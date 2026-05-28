cmake_minimum_required(VERSION 3.25)

project(IMASH LANGUAGES CXX)

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------
set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------
# Options
# ------------------------------------------------------------
# Enable this only when you want to run IMASH standalone tests
option(IMASH_BUILD_EXECUTABLE "Build IMASH test executable" ON)

# ------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)

# Look up IMAS-Cpp via pkg-config. PKG_CONFIG_PATH must already point at the
# directory containing al-cpp.pc (e.g. via an environment module, or set
# manually before invoking cmake).
pkg_check_modules(IMASCPP REQUIRED IMPORTED_TARGET al-cpp)

# ------------------------------------------------------------
# IMASH library
# ------------------------------------------------------------
add_library(IMASH STATIC
    src/IMASH.cpp
    src/TokamakData.cpp
)

target_include_directories(IMASH PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(IMASH PUBLIC
    PkgConfig::IMASCPP
    Threads::Threads
)

# Thread flags (needed by IMAS / al-cpp stack)
target_compile_options(IMASH PUBLIC -pthread)
target_link_options(IMASH PUBLIC -pthread)

# ------------------------------------------------------------
# Optional executable (for debugging / testing IMASH standalone)
# ------------------------------------------------------------
if(IMASH_BUILD_EXECUTABLE)

    add_executable(imash
        src/main.cpp
    )

    target_link_libraries(imash PRIVATE IMASH)

    # Ensure runtime can find IMAS shared libs
    set_target_properties(imash PROPERTIES
        BUILD_RPATH "${IMASCPP_LIBRARY_DIRS}"
        INSTALL_RPATH "${IMASCPP_LIBRARY_DIRS}"
        SKIP_BUILD_RPATH FALSE
    )

endif()

# ------------------------------------------------------------
# Debug / info output
# ------------------------------------------------------------
message(STATUS "---------------------------------------------")
message(STATUS "IMASCPP_LIBRARY_DIRS = ${IMASCPP_LIBRARY_DIRS}")
message(STATUS "IMASCPP_LIBRARIES    = ${IMASCPP_LIBRARIES}")
message(STATUS "IMASH_BUILD_EXECUTABLE = ${IMASH_BUILD_EXECUTABLE}")
message(STATUS "---------------------------------------------")
# ------------------------------------------------------------
# Install + CMake package export (added for EasyBuild packaging)
# ------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS IMASH
    EXPORT IMASHTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT IMASHTargets
    FILE IMASHTargets.cmake
    NAMESPACE IMASH::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IMASH
)

# Hand-written package config: re-resolves IMAS-Cpp via pkg-config on the consumer side,
# since the exported target depends on PkgConfig::IMASCPP which must be re-created.
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/IMASHConfig.cmake"
"include(CMakeFindDependencyMacro)
find_dependency(PkgConfig)
find_dependency(Threads)
pkg_check_modules(IMASCPP REQUIRED IMPORTED_TARGET al-cpp)
include(\"\${CMAKE_CURRENT_LIST_DIR}/IMASHTargets.cmake\")
")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/IMASHConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IMASH
)
