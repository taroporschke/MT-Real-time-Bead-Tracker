cmake_minimum_required(VERSION 3.10)
project(BeadTracker)

# Global settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenCV
set(OpenCV_DIR "C:/OpenCV/build")
find_package(OpenCV REQUIRED)

# DLL
add_library(BeadTrackerDLL SHARED
        src/BeadTrackerDLL.cpp
        src/BeadTrackerCore.cpp
        src/ForceCalculator.h
        src/test_main.cpp
        src/BeadTrackerCore.h
        src/BeadTrackerDLL.h
)

# Apply definitions and includes to the DLL
set_target_properties(BeadTrackerDLL PROPERTIES CXX_STANDARD 17)
target_compile_definitions(BeadTrackerDLL PRIVATE BEADTRACKERDLL_EXPORTS)
target_include_directories(BeadTrackerDLL PRIVATE src ${OpenCV_INCLUDE_DIRS})
target_link_libraries(BeadTrackerDLL ${OpenCV_LIBS})

# Text exe
add_executable(BeadTrackerTest
        src/test_main.cpp
)

# Include directories for the test
target_include_directories(BeadTrackerTest PRIVATE src ${OpenCV_INCLUDE_DIRS})

# Link executable to custom DLL and OpenCV
target_link_libraries(BeadTrackerTest BeadTrackerDLL ${OpenCV_LIBS})