add_executable(Testing "src/main.cpp")

target_include_directories(Testing PRIVATE "${PROJECT_SOURCE_DIR}/RedMath/src")

#target_link_libraries(Testing PRIVATE RedMath)