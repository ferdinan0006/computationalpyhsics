cmake_minimum_required(VERSION 3.20)

project(PenningTrap LANGUAGES CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Armadillo REQUIRED)

add_subdirectory(code)

find_package(Python3 COMPONENTS Interpreter REQUIRED)
set(VENV_DIR "${CMAKE_SOURCE_DIR}/.venv")

if(WIN32)
    set(VENV_PIP "${VENV_DIR}/Scripts/pip.exe")
else()
    set(VENV_PIP "${VENV_DIR}/bin/pip")
endif()

if(NOT EXISTS "${VENV_DIR}")
    message(STATUS "Creating Python venv at ${VENV_DIR}...")
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m venv "${VENV_DIR}"
        RESULT_VARIABLE VENV_RES
    )
    if(VENV_RES EQUAL 0 AND EXISTS "${CMAKE_SOURCE_DIR}/requirements.txt")
        message(STATUS "Installing requirements...")
        execute_process(
            COMMAND "${VENV_PIP}" install -r "${CMAKE_SOURCE_DIR}/requirements.txt"
        )
    endif()
endif()