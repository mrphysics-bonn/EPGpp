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

