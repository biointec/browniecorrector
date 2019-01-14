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
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile"
    exit
fi
if [ -z "$2" ]
  then
    echo "No argument supplied for coverage"
    echo "correct usage : ./louvain readFile  coverage iteration workDir staFile"
    exit
fi
if [ -z "$3" ]
  then
    echo "No argument supplied for iteration"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile"
    exit
fi
if [ -z "$4" ]
  then
    echo "No argument supplied for workDir"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile"
    exit
fi
if [ -z "$5" ]
  then
    echo "No argument supplied for staFile"
    echo "correct usage : ./louvain readFile coverage iteration workDir staFile"
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
##############

SOURCE="${BASH_SOURCE[0]}"
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
#############
i=0
networkFile=$workDir/network_0.dat
while [ $i -lt $iteration ]
do
   fileName="" 
   if [  -f $workDir/network_$i.dat  ]; then
        networkFile=$workDir/network_$i.dat
   fi
   $DIR/../release/src/convert -i $networkFile -o $workDir/graph.bin -w $workDir/graph.weights
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
   if [  -f $workDir/graph_node2comm_level$j ]; then
         rm $workDir/graph_node2comm_level$j  
   fi

   j=$[$j-1]
   k=1
   while [ $k -lt $j ]
   do
        rm $workDir/"graph_node2comm_level"$k
        k=$[$k+1]
   done
   if [  -f $workDir/graph.bin ]; then
         rm $workDir/graph.bin  
   fi
   if [  -f $workDir/graph.tree ]; then
         rm $workDir/graph.tree  
   fi
   if [  -f $workDir/graph.weights ]; then
         rm $workDir/graph.weights
   fi
   mv $fileName $workDir/$i.louvain
   gap=" "
   modularity=$(head -n 1 $workDir/result)
   if [  -f $workDir/result ]; then
         rm $workDir/result
   fi
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

     $DIR/../release/src/readclustering detectStableCommunity $workDir/"louvain.txt" $coverage $staFile $readNum $fullNames  
echo $DIR/../release/src/readclustering detectStableCommunity $workDir/"louvain.txt" $coverage $staFile $readNum $fullNames 
i=$[$i-1]
while [ $i -gt 0 ]
do
   if [  -f $workDir/network_$i.dat ]; then
        rm $workDir/network_$i.dat
   fi
   if [  -f $workDir/$i.louvain ]; then
        rm $workDir/$i.louvain
   fi
   i=$[$i-1]
done
fileName=$workDir/"louvain.txt"
echo $fileName
mkdir -p $workDir/clusters
rm -rf $workDir/clusters/*


