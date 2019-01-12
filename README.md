

# Short description:

Browniecorrector is a targeted error correction for short Illumina reads. unlike other error corrections which correct all the library, it focuses on the correction of only those reads that contain low complexity k-mers like homopolymers as such poly (A/ T). The main pipeline has four main steps:

1.selection of a highly repetitive kmer (default is a poly (A/T)).

1.extracts all reads that contain such a pattern 

2.clusters them into different groups using a community detection algorithm. 

3.consistant error correction in each cluster independently from other clusters. 

#  Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

#  Prerequisites

This package requires a number of packages to be install on your system. Required: CMake; Google's Sparsehash; gcc (GCC 4.7 or a more recent version) Optional: ZLIB; Googletest Unit Testing

How to install these packages:

As a root, execute the following commands:

on Redhat / Fedora distributions

    yum install cmake
    yum install sparsehash-devel
    yum install zlib-devel (optional)

on Ubuntu / Debian distributions

    aptitude install cmake
    aptitude install libsparsehash-dev
    aptitude install libghc-zlib-dev (optional)


# Installing browniecorrector:

The installation is now simple. First, clone browniecorrector from the Github address

    git clone "https://github.com/biointec/browniecorrector.git"

From this directory, run the following commands:

    mkdir release
    cd release
    cmake ..
    make install
The build directory should be named as the release, otherwise, you need to change that name in the pipeline scripts
By executing ./brownie you will see

    brownie: no command specified
    Try 'brownie --help' for more information

# Usage:

# Report bugs 
Please report bugs to : Mahdi.Heydari@UGent.be

# Citation

