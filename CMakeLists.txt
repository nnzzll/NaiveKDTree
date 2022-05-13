cmake_minimum_required(VERSION 3.0.0)
project(KDTree VERSION 0.1.0)

include(CTest)
enable_testing()

set(VTK_DIR "path/to/VTK-build")
find_package(VTK REQUIRED)
message("VTK VERSION:${VTK_VERSION}")
message("VTK_LIBRARIES:${VTK_LIBRARIES}")
message("VTK_INCLUDE_DIRS:${VTK_INCLUDE_DIRS}")
message("VTK_USE_FILE:${VTK_USE_FILE}")
include(${VTK_USE_FILE})

add_executable(KDTree main.cxx KDTree.h)
target_link_libraries(KDTree ${VTK_LIBRARIES})
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
