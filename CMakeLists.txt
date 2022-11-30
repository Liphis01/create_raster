cmake_minimum_required(VERSION 3.0.0)
project(main_exe VERSION 0.1.0)

include(CTest)
enable_testing()

file(
    GLOB_RECURSE
    source_files
    src/*.cpp
)

add_executable(
    ${PROJECT_NAME}
    ${source_files}
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

target_link_libraries(${PROJECT_NAME} proj)
include_directories("${PROJECT_SOURCE_DIR}/include")
