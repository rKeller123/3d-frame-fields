# CMakeList.txt : CMake project for 3d-odeco-frame-fields, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.10)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("3d-odeco-frame-fields")

# Add source to this project's executable.
add_executable (3d-odeco-frame-fields "src/3d-odeco-frame-fields.cpp" "src/3d-odeco-frame-fields.h"   "src/linalg.h"  "src/linalg.cpp" "src/odeco.h" "src/odeco.cpp")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET 3d-odeco-frame-fields PROPERTY CXX_STANDARD 20)
endif()


find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(3d-odeco-frame-fields PRIVATE Eigen3::Eigen)