add_subdirectory(lib)
add_subdirectory(tasks)

# Main executable
add_executable(${PROJECT_NAME}
    main.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE
    ${PROJECT_NAME}_tasks
    ${PROJECT_NAME}_lib
)