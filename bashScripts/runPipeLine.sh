#!/bin/bash


    #Copyright (C) <2019>  <mahdi.heydari@ugent.be>

    #This program is free software: you can redistribute it and/or modify
    #it under the terms of the GNU Affero General Public License as published
    #by the Free Software Foundation, either version 3 of the License, or
    #(at your option) any later version.

    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU Affero General Public License for more details.

    #You should have received a copy of the GNU Affero General Public License
    #along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    #

if [ -z "$1" ]
  then
    echo "No argument supplied for the input read file (fastq file) "
    echo "correct usage : ./brownieCorrection inputReadFile coverage tempDir "
    exit
fi
if [ -z "$2" ]
  then
    echo "No argument supplied for data coverage"
    echo "correct usage : ./brownieCorrection inputReadFile coverage tempDir "
    exit
fi
if [ -z "$3" ]
  then
    echo "No argument supplied for tempDir"
    echo "correct usage : ./brownieCorrection inputReadFile coverage tempDir "
    exit
fi
############ 0.1. Get input parameters
 inputLib=$1
 cov=$2
 workDir=$3 

############ 0.2. Set default parameters 

 iteration=20
 SOURCE="${BASH_SOURCE[0]}"
 DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
 mkdir -p $workDir 
 readFile=$workDir"/lowComplex.fastq"
 normalFile=$workDir"/normal.fastq"
 outputFile=$workDir"/brownie.fastq"
 KmerFile=$workDir"/kmerFile.fasta"
 
############ 1. kmer selection, kmer length should be and odd number. 
 pattern="AAAAAAAAAAAAAAA"
 #pattern="CCCCCCCCCCCCCCC"
 echo '>1' > $KmerFile
 echo $pattern >>$KmerFile
 kmer=${#pattern}
 echo $kmer
 if [ $((kmer%2)) -eq 0 ];
 then
    echo "The pattern length is used as kmer size and it cannot be an even number, change your pattern length";
    exit
 fi
 staFile=$workDir"/brownie."$kmer"/sta.txt" 

 ############ 2. Extract the reads that contain the given specific pattern. You need to provide the kmer file, the read file, kmer size and directory address.
 $DIR/../release/src/readclustering extractReads $KmerFile $inputLib $kmer $workDir
 
 ##### 3.1. Calclulate the similarity score between read pairs. 
 $DIR/calsimbrownieAligner.sh $readFile $kmer $cov $workDir $staFile
 
 ##### 3.2. Do the clustering with Louvain community detection algorithm. 

 $DIR/louvain.sh $readFile $cov $iteration $workDir $staFile
 
 fileName=$workDir/"louvain.txt"
 python3 $DIR/../pythonScripts/extractClusters.py $readFile $workDir/clusters $fileName $workDir/labelList.txt $workDir/network_0.dat
 
 ############  4.Correct each cluster with brownie
 $DIR/brownie.sh $workDir $kmer
 cat $normalFile $workDir/brownie.corrected.fastq >$outputFile
