cmake_minimum_required(VERSION 3.20)
project(QCL)

set(QRAM_LIB QRAMSimulator)

set(CMAKE_CXX_STANDARD 14)

#set(LLVM_DIR "E:\\llvm\\lib\\cmake\\llvm")
#set(MLIR_DIR "E:\\llvm\\lib\\cmake\\mlir")

find_package(MLIR REQUIRED CONFIG)

if (MLIR_FOUND)
    message(STATUS "Using MLIRConfig.cmake in: ${MLIR_DIR}")
    message(STATUS "Using LLVMConfig.cmake in ${LLVM_DIR}")

    set(LLVM_RUNTIME_OUTPUT_INTDIR ${CMAKE_BINARY_DIR}/bin)
    set(LLVM_LIBRARY_OUTPUT_INTDIR ${CMAKE_BINARY_DIR}/lib)

    list(APPEND CMAKE_MODULE_PATH "${MLIR_CMAKE_DIR}")
    list(APPEND CMAKE_MODULE_PATH "${LLVM_CMAKE_DIR}")

    include(TableGen)
    include(AddLLVM)
    include(AddMLIR)
    include(HandleLLVMOptions)

    include_directories(${LLVM_INCLUDE_DIRS})
    include_directories(${MLIR_INCLUDE_DIRS})

    link_directories(${LLVM_BUILD_LIBRARY_DIR})
    add_definitions(${LLVM_DEFINITIONS})
endif ()

include_directories(${PROJECT_SOURCE_DIR}/tpl
                    ${PROJECT_SOURCE_DIR}/tpl/argparse
                    ${PROJECT_SOURCE_DIR}/tpl/fmt/include)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

FIND_PACKAGE(OpenMP REQUIRED)

if (OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

add_subdirectory(frontend)
add_subdirectory(simulator)
add_subdirectory(tpl EXCLUDE_FROM_ALL)
add_executable(QCL main.cpp)
