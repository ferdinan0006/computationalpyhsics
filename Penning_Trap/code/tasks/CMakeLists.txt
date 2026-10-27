add_library(${PROJECT_NAME}_tasks
    src/Task8.cpp
    src/Task9.cpp
)

target_include_directories(${PROJECT_NAME}_tasks PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}_tasks PUBLIC
    ${PROJECT_NAME}_lib
)