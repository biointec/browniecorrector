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
if [  -f $workDir/brownie.corrected.fastq ]; then
        rm $workDir/brownie.corrected.fastq
fi
cd $workDir/clusters/
if [  -f $kmer.corrected.fastq ]; then
        rm $kmer.corrected.fastq
fi
kmerFile="kmer.txt"
SOURCE="${BASH_SOURCE[0]}"
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
newkmer=$kmer
trimming=1
echo "Wait until we correct all the clusters!!"
for i in {1..10000}
do
   if [  -f $i/$i.fastq  ]; then
     find ./$i/  ! -name $i.sam ! -name $i.fastq ! -name $i.sim ! -name karect_$i.fastq -type f -exec rm -rfd {} +
     kmerAnalysis=$i/"kmerAnalysis.txt"
     python3 $DIR/../pythonScripts/trimmer.py $i/$i.fastq $i/$i.trimmed.fastq >$kmerAnalysis
     cutoff=$(sed "1q;d" $kmerAnalysis )
     newkmer=$(sed "2q;d" $kmerAnalysis )
     outdir=$i/"mydir/"
     mkdir -p $outdir
     correctedFile=$outdir$i.cor.fastq 
    if [ "$trimming" -eq "1" ]; then
        /usr/bin/time -v $DIR/../release/src/brownie index $i/$i.trimmed.fastq -p $outdir -k $newkmer -nMM -nMEM >> $i/$i.brownieOut 2>>$i/$i.brownieOut_err
        if [  -f $outdir/arcs.stage4 ]; then
                rm $outdir/arcs.stage4
        fi
        /usr/bin/time -v $DIR/../release/src/brownie assemble -o $correctedFile $i/$i.fastq -p $outdir -k $newkmer -nMM -nMEM -c $cutoff >> $i/$i.brownieOut 2>>$i/$i.brownieOut_err
    else         
        /usr/bin/time -v $DIR/../release/src/brownie assemble -t 64 -p $outdir -k $newkmer -o $correctedFile $i/$i.fastq -nMM -nMEM -c $cutoff > $outdir/$i.brownieOut 2>$outdir/$i.brownieOut_err
    fi
    if [ -s "$correctedFile"  ]
    then
        cat  $correctedFile  >>$kmer.corrected.fastq
        #echo $correctedFile 
    else
        #echo $i":        unsuccessful error correction by brownie due to the low coverage of the reads in this cluster"
        cat  $i/$i.fastq  >>$kmer.corrected.fastq
     fi
   fi
done
cat  unclustered.fastq >>$kmer.corrected.fastq
mv $kmer.corrected.fastq ../brownie.corrected.fastq
cd $previousDir 


