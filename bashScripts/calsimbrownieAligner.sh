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
    echo "No argument supplied for the input low complexity read file (fastq file) "
    echo "correct usage : ./brownieCalSim readFile kmerSize coverage workDir staFile"
    exit
fi

if [ -z "$2" ]
  then
    echo "No argument supplied for kmer size (odd number )"
    echo "correct usage : ./brownieCalSim readFile kmerSize coverage workDir staFile"
    exit
fi
if [ -z "$3" ]
  then
    echo "No argument supplied for data coverage"
    echo "correct usage : ./brownieCalSim readFile kmerSize coverage tempDir staFile"
    exit
fi
if [ -z "$4" ]
  then
    echo "No argument supplied for workDir"
    echo "correct usage : ./brownieCalSim readFile kmerSize coverage workDir staFile"
    exit
fi
if [ -z "$5" ]
  then
    echo "No argument supplied for staFile"
    echo "correct usage : ./brownieCalSim readFile kmerSize coverage workDir staFile"
    exit
fi

   readFile=$1
   kmer=$2 
   extractKmer=$kmer
   cov=$3
   workDir=$4
   staFile=$5
   SOURCE="${BASH_SOURCE[0]}"
   DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
   
   baseFilename="lowComplex"
   lineNum=$(wc -l < $readFile)
   readNum=$((lineNum /8))
   mkdir -p $workDir/brownie.$kmer
   rm $workDir/brownie.$kmer/*
   
   $DIR/../release/src/brownie index -o $workDir/brownie.$kmer/$baseFilename.corr.fastq $readFile -p $workDir/brownie.$kmer/ -k $kmer -nMM -nMEM
   $DIR/../release/src/brownie align -o $workDir/brownie.$kmer/$baseFilename.corr.fastq $readFile -p $workDir/brownie.$kmer/ -k $kmer -nMM -nMEM
   rm $workDir/brownie.$kmer/$baseFilename.corr.fastq
   $DIR/../release/src/readclustering calSim $readFile  $workDir/brownie.$kmer/$baseFilename.ncf $workDir/network.dat $workDir/brownie.$kmer/nodes.stage5 $extractKmer $workDir/readNodeAss.txt $cov $staFile
   echo $DIR/../release/src/readclustering calSim $readFile $workDir/brownie.$kmer/$baseFilename.ncf $workDir/network.dat $workDir/brownie.$kmer/nodes.stage5 $extractKmer $workDir/readNodeAss.txt $cov $staFile
