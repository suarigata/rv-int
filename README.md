# RV-INT: The Risc-V Interpreter

RV-INT implements a fast Risc-V interpreter.

## Bulding it

To build RV-INT, you are going to need g++7.1 or later with C++17 support, CMake 2.8 or later, PAPI and LLVM 6.0. After installing all dependencies, it is a simple cmake/make usage:

```
cd rv-int
mkdir build
cd build
cmake ..
make
sudo make install
```

### Usage

After building and installing rv-int, you can easily use it to emulate Risc-V elf binaries using the following commands:

```
rv-int -bin PathToBinary [-v] [-args <args>]

ARGUMENTS:
  -bin : Path to the binary which will be emulated.
  -h : Displays the help message
  -v : Displays the OpenISA instructions from the compiled regions
```

### LICENSE

This project is being developed at the Institute of Computing - Unicamp as part of @suarigata undergraduate final project. You are free to contact him and use this code under the MIT LICENSE.
