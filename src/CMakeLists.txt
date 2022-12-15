file(
    GLOB_RECURSE
    source_files
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_library(
    ${PROJECT_NAME}
    SHARED
    ${source_files}
)

target_link_libraries(${PROJECT_NAME} proj)
include_directories("${PROJECT_SOURCE_DIR}/include")