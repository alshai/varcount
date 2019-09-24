cmake_minimum_required(VERSION 3.10)

project(varcount CXX C)

set(PROJECT_URL "https://github.com/alshai/varcount")
set(PROJECT_DESCRIPTION "Outputs alignment coverage over a predefined set of variants")
set(CMAKE_CXX_STANDARD 11)

# htslib
find_package(PkgConfig REQUIRED)
pkg_search_module(PKG_HTSLIB REQUIRED htslib)

add_executable(varcount varcount.cpp hts_util.hpp)
target_link_libraries(varcount ${PKG_HTSLIB_LDFLAGS})

# add_executable(next_gt next_gt.cpp hts_util.hpp)
# target_link_libraries(next_gt ${PKG_HTSLIB_LDFLAGS})

INSTALL(TARGETS varcount DESTINATION bin)
