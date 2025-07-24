cmake_minimum_required(VERSION 3.21)

project(task1)
set(cmake_cxx_standart_required 23)

find_package(GTest REQUIRED)

# add_compile_options(-fsanitize=address,undefined,leak -fno-sanitize-recover=all)
# add_link_options(-fsanitize=address,undefined,leak)

file(GLOB HEADERS src/*.h)
file(GLOB SOURCES src/*.cpp)

file(GLOB TEST_SOURCES tests/*.cpp)

add_executable(task-exec main.cpp ${SOURCES})
add_executable(tests ${TEST_SOURCES} ${HEADERS} ${SOURCES})

target_include_directories(tests PRIVATE src tests)
target_link_libraries(tests PUBLIC GTest::gtest GTest::gtest_main)
target_include_directories(task-exec PRIVATE "${PROJECT_SOURCE_DIR}/src")