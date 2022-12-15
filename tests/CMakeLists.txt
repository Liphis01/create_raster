include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  FetchContent_MakeAvailable(googletest)
  
add_executable(test_shadows test_shadows.cpp)

target_link_libraries(test_shadows GTest::gtest_main)
target_link_libraries(test_shadows ${PROJECT_NAME})

include_directories("${PROJECT_SOURCE_DIR}/include")

include(GoogleTest)
gtest_discover_tests(test_shadows)

add_test(
    NAME test_shadows
    COMMAND test_shadows
)