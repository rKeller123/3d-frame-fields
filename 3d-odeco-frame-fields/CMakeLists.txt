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

add_library(core OBJECT
    src/linalg.cpp
    src/odeco.cpp
    src/newton.cpp
)

target_include_directories(core PUBLIC src)
target_link_libraries(core PUBLIC Eigen3::Eigen)
set_target_properties(core PROPERTIES CXX_STANDARD 20)

target_compile_options(core PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:
        -O3
        -march=native
        -DEIGEN_NO_DEBUG
    >
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:
        /O2
        /arch:AVX2
        /DEIGEN_NO_DEBUG
    >
)

# Add source to this project's executable.
add_executable(3d-odeco-frame-fields src/3d-odeco-frame-fields.cpp)
add_executable(cli src/cli.cpp)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET 3d-odeco-frame-fields PROPERTY CXX_STANDARD 20)
  set_property(TARGET cli PROPERTY CXX_STANDARD 20)
endif()

target_compile_options(cli PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:
        -O3
        -march=native
        -DEIGEN_NO_DEBUG
    >
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:
        /O2
        /arch:AVX2
        /DEIGEN_NO_DEBUG
    >
)
target_compile_options(3d-odeco-frame-fields PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:
        -O3
        -march=native
        -DEIGEN_NO_DEBUG
    >
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:
        /O2
        /arch:AVX2
        /DEIGEN_NO_DEBUG
    >
)


find_package(Eigen3 CONFIG REQUIRED)
target_link_libraries(3d-odeco-frame-fields PRIVATE core Eigen3::Eigen)
target_link_libraries(cli PRIVATE core Eigen3::Eigen)