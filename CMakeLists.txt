cmake_minimum_required(VERSION 3.20)
project(NumLab LANGUAGES CXX)

# --- ustawienia kompilatora ---------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Gdy u¿ytkownik nie poda -DCMAKE_BUILD_TYPE, wymuœ Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# --- Ÿród³a biblioteki ---------------------------------------------------
file(GLOB LIB_SOURCES "src/*.cpp")

add_library(NumLab STATIC ${LIB_SOURCES})

target_include_directories(NumLab PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# --- instalacja (opcjonalna) --------------------------------------------
install(TARGETS NumLab ARCHIVE DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
