

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
To run BrownieCorrector, you need to run the runPipeLine.sh script in the bashScrips folder. You need to provide the script three arguments like this :

    ./runPipeLine.sh  inputReadFile coverage tempDir
    
The first argument is inputReadFile is the Fastq file library in which the read pairs are interleaved.  Therefore two consecutive reads (one pair) must have an identical id. In case they have different tags, for example, @firstID/1 and @secondID/2,  remove the /1 and /2.  Also, please make sure that the file contains reads from a single library (same insert size). 

The second argument is the coverage of the library. If you don't have the coverage but you know the approximate genome size you can compute the coverage like this: C = (L*N)/G. C stands for coverage, G is the genome length, L is the read length, and N is the number of reads. 


The third argument is the working directory to keep the temporary files and the results.

Inside the runPipeLine.sh file you can change some parameters. For example, you can change the pattern from "AAAAAAAAAAAAAAA" to other patterns. Please make sure that the length of this pattern should be an odd number. The pipeline has four states :

1. Extract reads:  it loops over the library and extracts all the pairs that contain the given pattern. The input of this step is the read library (first argument given to run the pipeline) and the output is two files: lowComplex.fastq and normal.fastq. BrownieCorrector doesn't touch the normal file, it only corrects the reads in lowComplex.fastq. 
2. Compute similarity score: it computes the similarity score between all pairs inside lowComplex.fastq file. Overlap alignment is used to find the similarity score. Because this process is quadratic and time-consuming (finding the similarity score between all pairs), we only compute the similarity score if they share a mutual k-mer. We first align the reads to the graph, two pair that align to the same node share a mutual k-mer. But of course, all the pairs share at least one node (node that contains the pattern). Based on the code coverage, we only take into account nodes whose multiplicity is 1. This section can be replaced by your method of choice. The output of this step is a tab-delimited file (network_0.dat) that shows the similarity score between every two pairs. Those that have a similarity score lower than some threshold are filtered. Therefore, the file contains a line like this:

     
     First Pair  secodn Pair score
     1	24236	32.00
     1	39078	70.50
     1	70713	33.50
     1	70879	27.50
     1	73169	32.50
     1	77762	23.00

Reads are numbered based on their order in lowComplex.fastq file. So 1 point to the two first reads in that file (one pair). You can use a different way to calculate the similarity score, but remember your output file should be similar to this file.

3. Read Clustering: it finds the clusters based on the similarity between pairs. We used the Louvain community detection algorithm to do the clustering. The input of this step is the output of the previous step (network_0.dat), and the output is the tab-deliminate file (louvain.txt) that shows each pair belongs to which cluster. For example, there are few lines of this file: (for example, the first and the third read are in the same cluster)

    read-Pair   clusterID
    1	1
    2	1007
    3	1
    4	1319

Based on this clustering, reads are assigned to different groups and different folders are made for each group in clusters folder. 

4.Error correction. In this step, reads in the cluster are corrected independently. The corrected reads in each cluster are concatenated to each other in brownie.corrected.fastq file. Then brownie.corrected.fastq and normal are concatenated together and the result is brownie.fastq which is the output of the correction. Remember, BrownieCorrector doesn't keep the order of reads as given. 



# Report bugs 
Please report bugs to : Mahdi.Heydari@UGent.be

# Citation

