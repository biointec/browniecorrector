

#GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
#Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.


import sys
from Bio import SeqIO
from statistics import mean, median , pstdev
import os
import shutil
minClusterLen = 10
maxClusterLen = 300

def RC(seq):
    seq_dict={'A':'T', 'T':'A', 'G':'C','C':'G', 'N':'N','a':'t', 't':'a', 'g':'c','c':'g', '-':'-'}
    return "".join([seq_dict[base] for base in reversed(seq)])
def createCluster(directory, fileName):
    directory = directory + "/"+ str (fileName)
    if not os.path.exists(directory):
        os.makedirs(directory)
        #print ("we made this directory" + directory)
    directory = directory +  "/" + str (fileName) + ".fastq"
    f = open (directory, 'w')
    f.close()
def extractClusterInfo(comunitiesFileName, directory):
    f = open (comunitiesFileName)
    content  =f.read()
    items = content.split("\n")
    items = items[:-1]
    itemsInfo = {}
    numOfItems = len (items)
    freqList = {}
    for item in items:
        elements = item.split()
        try :        
            freqList[int (elements[1])] = 1 +  freqList[int (elements[1])]
        except :
            freqList[int (elements[1])] = 1
            
    clusters =set()
    clusterSize=[]
    readsinHigherSizeCluster=0
    readsinLowerSizeCluster=0
    for item in items:
        elements = item.split()
        elementSet = set ()
        if (freqList[int (elements[1])]  >=minClusterLen  and freqList[int (elements[1])] < maxClusterLen):
            clusters.add(int (elements[1]))
            clusterSize.append(freqList[int (elements[1])])
            createCluster(directory, int (elements[1]))
            
            for element in elements :
                itemsInfo[int (elements[0])] = int (elements[1])
        else:
            if (freqList[int (elements[1])] > maxClusterLen):
                readsinHigherSizeCluster  = readsinHigherSizeCluster +1
            if (freqList[int (elements[1])] < minClusterLen):
                readsinLowerSizeCluster = readsinLowerSizeCluster+1
                
    print ("Number of reads in clusters lower than threshold:\t"+ str (readsinLowerSizeCluster))
    print ("Number of reads in clusters higher than threshold:\t"+ str (readsinHigherSizeCluster))
    print ("Total number of unclustered reads:\t"+str (readsinLowerSizeCluster+readsinHigherSizeCluster) )
    print ("Number of clusters : " +str ( len (clusters) ))
    avg_value = mean(clusterSize)
    median_value = median(clusterSize)
    print ("Average size of clusters: ",round( avg_value))
    print ("Median size of clusters" , median_value)
    
    return itemsInfo
def writeReadsIntoClustersWithInternalSim(inputFileName, comunitiesFileName, directory, labelListFileName, networkFileName):
    dic = extractClusterInfo(comunitiesFileName, directory)
    labelList = open(labelListFileName, "w")
    clusters={}
    restFile = open (directory + "/unclustered.fastq" , "w")
    labelList = open(labelListFileName, 'w')
    pairNum = 1
    num = 1
    firstID = ""
    firstRead = ""
    firstQuality =""
    namesIndex ={}
    with open(inputFileName, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):

            if (num%2==1):
                firstID = record.id
                firstRead = str (record.seq)
                firstQuality = str(record.format("fastq")).split("\n")[3]
                
            else:
                secondID = record.id
                secondRead = str (record.seq)
                secondQuality = str(record.format("fastq")).split("\n")[3]
                try :
                    namesIndex [pairNum] = record.id.split("/")[0]
                    clusterID = dic [pairNum]
                    try :
                        clusters[clusterID].add(pairNum)
                    except :
                        readSet = set()
                        readSet.add(pairNum)
                        clusters[clusterID] = readSet
                        
                    fileName = directory + "/"+  str (clusterID)+"/"+ str (clusterID) +".fastq"
                    f = open (fileName, 'a')
                    f.write("@"+str(firstID)+"\n")
                    f.write(str(firstRead)+"\n")
                    f.write("+\n")
                    f.write(str(firstQuality)+"\n")
                    f.write("@"+str(secondID)+"\n")
                    f.write(secondRead+"\n")
                    f.write("+\n")
                    f.write( secondQuality+"\n")

                    
                except :
                    clusterID=0
                    restFile.write("@"+str(firstID)+"\n")
                    restFile.write(str(firstRead)+"\n")
                    restFile.write("+\n")
                    restFile.write(str(firstQuality)+"\n")
                    restFile.write("@"+str(secondID)+"\n")
                    restFile.write(secondRead+"\n")
                    restFile.write("+\n")
                    restFile.write( secondQuality+"\n")
                labelList.write(str (pairNum)+"\t"+str(clusterID)+"\n")    

                pairNum = pairNum+ 1
                
                
            num = num +1
    num = (num - 1)/2+1
    labelList.close()
    print (networkFileName)
    wirteInternalSim(networkFileName, clusters, num, namesIndex, directory)
def writeReadsIntoClusters(inputFileName, comunitiesFileName, directory, labelListFileName):
    labelList = open(labelListFileName, "w")
    dic = extractClusterInfo(comunitiesFileName, directory)
    clusters={}
    restFile = open (directory + "/unclustered.fastq" , "w")
    pairNum = 1
    num = 1
    firstID = ""
    firstRead = ""
    firstQuality =""
    namesIndex ={}
    overcorrectedClusters = {}
    with open(inputFileName, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):

            if (num%2==1):
                firstID = record.id
                firstRead = str (record.seq)
                firstQuality = str(record.format("fastq")).split("\n")[3]
                firstIDFullLine = str(record.format("fastq")).split("\n")[0]
            else:
                secondID = record.id
                secondRead = str (record.seq)
                secondQuality = str(record.format("fastq")).split("\n")[3]
                secondIDFullLine = str(record.format("fastq")).split("\n")[0]

                try :
                    namesIndex [pairNum] = record.id.split("/")[0]
                    clusterID = dic [pairNum]
                    try :
                        clusters[clusterID].add(pairNum)
                    except :
                        readSet = set()
                        readSet.add(pairNum)
                        clusters[clusterID] = readSet
                    if ( firstIDFullLine.find("A(1)A")>-1 or secondIDFullLine.find("A(1)A")>-1):
                        try:
                            overcorrectedClusters[clusterID] = overcorrectedClusters[clusterID]+1
                        except:
                            overcorrectedClusters[clusterID] = 1
                    
                    fileName = directory + "/"+  str (clusterID)+"/"+ str (clusterID) +".fastq"
                    f = open (fileName, 'a')
                    f.write(firstIDFullLine+"\n")
                    f.write(str(firstRead)+"\n")
                    f.write("+\n")
                    f.write(str(firstQuality)+"\n")
                    f.write(secondIDFullLine+"\n")
                    f.write(secondRead+"\n")
                    f.write("+\n")
                    f.write( secondQuality+"\n")
                    f.close()
                except :
                    clusterID=0
                    restFile.write(firstIDFullLine+"\n")
                    restFile.write(str(firstRead)+"\n")
                    restFile.write("+\n")
                    restFile.write(str(firstQuality)+"\n")
                    restFile.write(secondIDFullLine+"\n")
                    restFile.write(secondRead+"\n")
                    restFile.write("+\n")
                    restFile.write( secondQuality+"\n")
                labelList.write(str (pairNum)+"\t"+str(clusterID)+"\n")    
                pairNum = pairNum+ 1
            num = num +1
    labelList.close();
    num = (num - 1)/2+1
    for key, value in overcorrectedClusters.items():
        fileName = directory + "/"+  str (key)+"/"+ str (key) +".fastq"
        with open(fileName, "rU") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                    restFile.write("@"+str(record.id)+"\n")
                    restFile.write(str(record.seq)+"\n")
                    restFile.write("+\n")
                    restFile.write(str(str(record.format("fastq")).split("\n")[3])+"\n")
        shutil.rmtree(directory + "/"+  str (key)+"/")
        
    
def wirteInternalSim(networkFileName, clusters, numOfPairs, namesIndex, directory):

    networkFile = open(networkFileName)
    
    Matrix = {}
    for line in networkFile:
        elements = line.split()
        source =int ( elements[0] )
        target =int ( elements[1] )
        sim = float (elements[2])
        if (source <= target):
            Matrix[source, target] = sim
        else:
            Matrix[target, source] = sim

    for clusterID, readSet in clusters.items():
        fileName = directory + "/"+ str (clusterID) + "/" + str (clusterID)  +".sim"
        f = open (fileName, 'w')
        f.write("Source \t Target \t Weight \t extra \n")
        readList = list (readSet)
        i = 0
        while (i < len (readList)):
            j = i;
            source =readList[i]
            sourceName = namesIndex[source]
            while (j< len (readList)):
                target =readList[j]
                targetName = namesIndex[target]
                try :
                    sim = 0
                    if (source <=target):
                        sim =  Matrix [source, target]
                    else:
                        sim =  Matrix [target, source]
                    if (sim>0):    
                        f.write (str (source) + "\t" +str (target) + "\t" +str ( sim)+  "\t(" +str (sourceName) + ","+ str (targetName) +")\n")
                    j = j+1
                except :
                    #print (source , target, 0)
                    j = j+1
                    continue;
            i = i +1
        f.close()

    
def main(argv=None):
    if argv == None:
        argv = sys.argv
    if len(argv) < 1:
        exit()
    readFile = str(sys.argv[1])
    print ("read file :"+ readFile)
    directory = str(sys.argv[2])
    print ("directory :"+ directory)
    comunitiesFileName = str(sys.argv[3])
    print ("comunitiesFileName :"+ comunitiesFileName)
    labelListFileName=str(sys.argv[4])
    print ("labelListFileName :"+ labelListFileName)
    if (len(argv) == 6):
        print ("network sim file is given")
	
        networkFile =  str(sys.argv[5])
        writeReadsIntoClustersWithInternalSim(readFile, comunitiesFileName , directory,labelListFileName, networkFile)
    else:
        print ("network sim file is not given")
        writeReadsIntoClusters(readFile, comunitiesFileName , directory,labelListFileName)

main() 






