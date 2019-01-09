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
    echo "correct usage : ./brownie workDir kmerSize "
    exit
fi
if [ -z "$2" ]
  then
    echo "No argument supplied for kmer size"
    echo "correct usage : ./brownie workDir kmerSize  "
    exit
fi
workDir=$1
kmer=$2  # the second argument is the kmer size

previousDir=$(pwd)

rm $workDir/brownie.corrected.fastq
cd $workDir/clusters/
rm $kmer.corrected.fastq
kmerFile="kmer.txt"
SOURCE="${BASH_SOURCE[0]}"
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
newkmer=$kmer
echo "Wait until we correct all the clusters!!"
for i in {1..10000}
do
   if [  -f $i/$i.fastq  ]; then
     find ./$i/  ! -name $i.sam ! -name $i.fastq ! -name $i.sim ! -name karect_$i.fastq -type f -exec rm -rfd {} +
    kmerAnalysis=$i/"kmerAnalysis.txt"
    correctedFileTrimmed=$i/trimmed/$i.cor.fastq
    python3 $DIR/../pythonScripts/trimmer.py $i/$i.fastq $i/$i.trimmed.fastq >$kmerAnalysis
    cutoff=$(sed "1q;d" $kmerAnalysis )
    mkdir -p  $i/trimmed
    /usr/bin/time -v $DIR/../release/src/brownie index $i/$i.trimmed.fastq -p $i/trimmed -k $kmer -nMM -nMEM >> $i/$i.brownieOut 2>>$i/$i.brownieOut_err
    rm $i/trimmed/arcs.stage4
    /usr/bin/time -v $DIR/../release/src/brownie assemble -o $correctedFileTrimmed $i/$i.fastq -p $i/trimmed/ -k $kmer -nMM -nMEM -c $cutoff >> $i/$i.brownieOut 2>>$i/$i.brownieOut_err
    if [ -s "$correctedFileTrimmed"  ]
    then
        cat  $correctedFileTrimmed  >>$kmer.corrected.fastq
    else
        #echo $i":        unsuccessful error correction by brownie due to the low coverage of the reads in this cluster"
        cat  $i/$i.fastq  >>$kmer.corrected.fastq
     fi
   fi
done
cat  unclustered.fastq >>$kmer.corrected.fastq
mv $kmer.corrected.fastq ../brownie.corrected.fastq
cd $previousDir 

