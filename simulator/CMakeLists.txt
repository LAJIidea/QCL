file(GLOB_RECURSE SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_library(${QRAM_LIB} STATIC ${SRC})

target_include_directories(${QRAM_LIB} PUBLIC  .)