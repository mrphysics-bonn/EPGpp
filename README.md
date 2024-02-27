# EPGpp: C++ implementation of Extended Phase Graphs with Cython wrapper 
Authors: Tony Stöcker, Ed Pracht, DZNE Bonn


## Installation 
### C++ part 
(assuming you are in the projects main directory)
```
    mkdir build
    cd build
    cmake ..
    make
    make install
    cd ..
```
If write permission to /usr/local is not available, install directory can be defined in the command line. For example:
```
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/any_defined_directory ..
    make
    make install
    cd ..
```
and comment out ```SET(CMAKE_INSTALL_PREFIX $ENV{HOME}/.local)``` in the CMakeLists.txt in main directory.
The install is optional and not needed for the Cython wrapper

### Cython wrapper 
Get Cython first:
```
    conda install cython
```
now, assuming you are in the projects main directory:
```
    pip install .
```
## Usage
Have a look at the pure C++ example (source: main.cpp) first:
```
    cd build
    ./epg
```
Afterwards run the Jupyter notebook (example.ipynb) to see how the Cython wrapper works.

