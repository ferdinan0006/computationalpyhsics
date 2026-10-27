add_library(${PROJECT_NAME}_lib
    src/Particle.cpp
    src/PenningTrap.cpp
    src/PenningTrapPerturbed.cpp
    src/Solvers.cpp
)

target_include_directories(${PROJECT_NAME}_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ARMADILLO_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME}_lib PUBLIC
    ${ARMADILLO_LIBRARIES}
)