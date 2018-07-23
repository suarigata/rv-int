cmake_minimum_required(VERSION 3.5)
project(riscv-int)

set(CMAKE_CXX_FLAGS "-g -Ofast -Werror -Wall -Wextra -Wno-unused-function -std=c++17 -fno-rtti -DNDEBUG")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Ofast")

include_directories(include)

add_subdirectory(arglib)

add_library(int interpreter.cpp machine.cpp syscall.cpp RVDecoder.cpp)

add_executable(oi-int main.cpp)

install(TARGETS oi-int RUNTIME DESTINATION bin)

target_link_libraries(oi-int
  int 
  arglib 
  z
)
