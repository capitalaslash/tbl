cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(tbl
  VERSION 0.1
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic")

include(FetchContent)
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

include(FeatureSummary)

# eigen
find_package(Eigen3 REQUIRED NO_MODULE)

# libmesh
find_package(PkgConfig REQUIRED)
pkg_check_modules(libmesh REQUIRED IMPORTED_TARGET libmesh)
# workaround for bogus libmesh include path
get_target_property(
  LIBMESH_INCLUDE_DIRECTORIES PkgConfig::libmesh INTERFACE_INCLUDE_DIRECTORIES)
list(REMOVE_ITEM LIBMESH_INCLUDE_DIRECTORIES "/usr/include/tirpc")
set_property(TARGET PkgConfig::libmesh
  PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${LIBMESH_INCLUDE_DIRECTORIES})

# igl
find_package(libigl REQUIRED)

# cgal
list(APPEND CMAKE_MODULE_PATH "$ENV{SPACK_ENV}/.spack-env/view/lib/cmake/CGAL")
list(APPEND CMAKE_MODULE_PATH "$ENV{SPACK_ENV}/.spack-env/view/lib64/cmake/CGAL")
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
include(CGALConfig)
include(UseCGAL)

add_executable(main)

set(SRC_FILES
  main.cpp
)
# set(HEADER_FILES
#   header.hpp
# )

target_sources(main
	PRIVATE ${SRC_FILES}
	# PUBLIC ${HEADER_FILES}
)
target_link_libraries(main
  PUBLIC
    Eigen3::Eigen
    PkgConfig::libmesh
    igl::igl_core
    CGAL::CGAL
    fmt::fmt
)

feature_summary(WHAT ALL)

find_program(CLANGFORMAT clang-format)
if (CLANGFORMAT)
  add_custom_target(format
    COMMAND ${CLANGFORMAT}
    --verbose -i
    ${SRC_FILES} ${HEADER_FILES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )
else()
  add_custom_target(format
    COMMAND echo "clang-format could not be found"
  )
endif()

