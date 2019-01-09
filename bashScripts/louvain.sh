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
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile clusteringfileName"
    exit
fi
if [ -z "$2" ]
  then
    echo "No argument supplied for coverage"
    echo "correct usage : ./louvain readFile  coverage iteration workDir staFile clusteringfileName"
    exit
fi
if [ -z "$3" ]
  then
    echo "No argument supplied for iteration"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile clusteringfileName"
    exit
fi
if [ -z "$4" ]
  then
    echo "No argument supplied for workDir"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile clusteringfileName"
    exit
fi
if [ -z "$5" ]
  then
    echo "No argument supplied for staFile"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile clusteringfileName"
    exit
fi
if [ -z "$6" ]
  then
    echo "No argument supplied for clusteringfileName"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile clusteringfileName"
    exit
fi


#################
readFile=$1
baseFilename="${readFile%.*}"
fullNames=""
coverage=$2
iteration=$3
workDir=$4
staFile=$5
clusteringfileName=$6
##############

SOURCE="${BASH_SOURCE[0]}"
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
#############
i=0
while [ $i -lt $iteration ]
do
   fileName="" 
   $DIR/../release/src/convert -i $workDir/network_$i.dat -o $workDir/graph.bin -w $workDir/graph.weights
   $DIR/../release/src/louvain   $workDir/graph.bin -l -1 -q 0 -w $workDir/graph.weights > $workDir/graph.tree 2>$workDir/result
   j=1
   $DIR/../release/src/hierarchy $workDir/graph.tree -l $j > $workDir/graph_node2comm_level$j 
   while [ -s $workDir/graph_node2comm_level$j ]
   do
        fileName=$workDir/graph_node2comm_level$j 
        sed -i -e "1d" $fileName
        echo $fileName
        j=$[$j+1]
        $DIR/../release/src/hierarchy $workDir/graph.tree -l $j > $workDir/graph_node2comm_level$j 
   done
   rm $workDir/graph_node2comm_level$j
   j=$[$j-1]
   k=1
   while [ $k -lt $j ]
   do
        rm $workDir/"graph_node2comm_level"$k
        k=$[$k+1]
   done
   rm $workDir/graph.bin
   rm $workDir/graph.tree
   rm $workDir/graph.weights
   mv $fileName $workDir/$i.louvain
   gap=" "
   modularity=$(head -n 1 $workDir/result)
   rm $workDir/result
   minModularity=.5
   st=`echo "$modularity > $minModularity" | bc`
   if [ $st -eq 1 ] ; then
       fullNames=$fullNames$workDir/$i.louvain$gap
       echo $modularity
   fi
   i=$[$i+1]
done
lineNum=$(wc -l < $readFile)
readNum=$((lineNum /8))
konsule_output=$workDir/"connectedComp_output"
konsule_output_err=$workDir/"connectedComp_err"

     $DIR/../release/src/readclustering detectStableCommunity $clusteringfileName $coverage $staFile $readNum $fullNames  
echo $DIR/../release/src/readclustering detectStableCommunity $clusteringfileName $coverage $staFile $readNum $fullNames 
i=1
while [ $i -lt $iteration ]
do
   rm $workDir/network_$i.dat
   rm $workDir/$i.louvain
   i=$[$i+1]
done

mkdir -p $workDir/clusters
rm -rf $workDir/clusters/*


