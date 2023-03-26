# MP4GNN 

MP4GNN is an MPI-based library for partitioning large irregular graphs. MP4GNN is based on [ParMETIS'](https://github.com/KarypisLab/ParMETIS) multilevel multi-constraint k-way 
algorithm developed in our lab and contains various efficiency and memory
optimizations. Many of these optimizations were developed while working on
[DistDGL](https://arxiv.org/abs/2010.05337) and
[DistDGLv2](https://arxiv.org/abs/2112.15345).

##  Downloading MP4GNN 

You can download MP4GNN by simply cloning it using the command:
```
git clone git@github.com:KarypisLab/PM4GNN.git
```

## Building the MP4GNN library

To build MP4GNN you can follow the instructions below:

### Dependencies

General dependencies for building MP4GNN are: gcc, cmake, build-essential, and an MPI library. 
In Ubuntu systems these can be obtained from the apt package manager (e.g., apt-get install cmake, mpich, etc) 

```
sudo apt-get install build-essential
sudo apt-get install cmake
```

In addition, you need to download and install
[GKlib](https://github.com/KarypisLab/GKlib) and 
[METIS](https://github.com/KarypisLab/METIS) by following the instructions there. 


### Building and installing MP4GNN  

MP4GNN is primarily configured by passing options to make config. For example:

```
make config cc=mpicc prefix=~/local
make install
```

will configure MP4GNN to be built using mpicc and then install the binaries, header files, and libraries at 

```
~/local/bin
~/local/include
~/local/lib
```

directories, respectively.

### Common configuration options are:

    cc=[compiler]     - The C compiler to use [default is determined by CMake]
    shared=1          - Build a shared library instead of a static one [off by default]
    prefix=[PATH]     - Set the installation prefix [~/local by default]
    gklib_path=[PATH] - Set the prefix path where GKlib has been installed. You can skip
                        this if GKlib's installation prefix is the same as that of
                        MP4GNN.
    metis_path=[PATH] - Set the prefix path where METIS has been installed. You can skip
                        this if METIS' installation prefix is the same as that of
                        MP4GNN.

### Advanced debugging related options:

    gdb=1           - Build with support for GDB [off by default]
    debug=1         - Enable debugging support [off by default]
    assert=1        - Enable asserts [off by default]
    assert2=1       - Enable very expensive asserts [off by default]

### Other make commands

    make uninstall
         Removes all files installed by 'make install'.

    make clean
         Removes all object files but retains the configuration options.

    make distclean
         Performs clean and completely removes the build directory.


### Definitions of supported data types

MP4GNN uses the same data types for integers and floating point numbers (32/64 bit
integers and single/double precision floating point numbers) as used when configuring
and building METIS.


## Copyright & License Notice
Copyright 1998-2023, Regents of the University of Minnesota and is licensed under the 
Apache License, Version 2.0. 


## Cite 

There are several papers describing ParMetis' underlying algorithms
[[1]](https://dl.acm.org/doi/abs/10.1145/369028.369103)
[[2]](https://www.sciencedirect.com/science/article/abs/pii/S0743731597914040)
[[3]](https://www.sciencedirect.com/science/article/abs/pii/S0743731597914039)
[[4]](https://link.springer.com/chapter/10.1007/3-540-44520-X_39)

For citing this particular repository, use the following:

```
@misc{pm4gnn23,
  author = {George Karypis},
  title = {PM4GNN -- Parallel Graph Partitioning for Distributed GNNs},
  howpublished = "\url{https://github.com/KarypisLab/PM4GNN/}",
  year = {2023}
}
```
