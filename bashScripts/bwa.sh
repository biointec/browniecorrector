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
    echo "correct usage : ./bwa workDir genome.fasta  "
    exit
fi

if [ -z "$2" ]
  then
    echo "No argument supplied for genome file (fasta file)"
    echo "correct usage : ./bwa workDir genome.fasta "
    exit
fi


workDir=$1
genomeFile=$2

fileName=$workDir/"louvain.txt"
genomeAdd=$(readlink -f $genomeFile)

previousDir=$(pwd)
cd $workDir/clusters/
for i in {1..1000000}
do

    if [  -f $i/$i.fastq  ]; then
    cd $i
    bwa mem $genomeAdd $i.fastq -p -t 4 | samtools view  -Sb - | samtools sort -o sorted.bam && samtools index sorted.bam | samtools view -h -o $i.sam sorted.bam    
    #bwa mem $genomeAdd $i.fastq -p -t 4 | samtools view  -Sb - | samtools sort - sorted && samtools index sorted.bam | samtools view -h -o $i.sam sorted.bam 
    rm sorted.bam
    rm sorted.bam.bai
    cd ..
   fi
done

cd $previousDir
