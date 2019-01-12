    #Copyright (C) <2019> Mahdi Heydari <mahdi.heydari@ugent.be>

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
import sys
from Bio import SeqIO
from statistics import mean, median , pstdev
import os
import shutil
import math
kmerSize = 41
minCutOff=1
maxCutoff=5
minKmerSize = 15
maxKmerSize = 41
minCov=5
minGenomeSize=100
def RC(seq):
    seq_dict={'A':'T', 'T':'A', 'G':'C','C':'G', 'N':'N','a':'t', 't':'a', 'g':'c','c':'g', '-':'-'}
    return "".join([seq_dict[base] for base in reversed(seq)])
def trimReads(inputFileName , outputFileName):
    outfile =open (outputFileName, 'w')
    with open(inputFileName, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            firstID = record.id
            firstRead = str (record.seq)
            firstQuality = str(record.format("fastq")).split("\n")[3]
            i = len (firstQuality)-1
            while (i>= 0 and firstQuality[i]<='#'):
                i = i -1
            newRead = firstRead[0:i+1]
            newQual=firstQuality[0:i+1]
            if (len (newRead) > 0):
                outfile.write ("@"+firstID+"\n")
                outfile.write(newRead+"\n")
                outfile.write("+\n")
                outfile.write(newQual+"\n")
    outfile.close()
def computeAvg(GToneFreq):
    sum = 0.0
    num = 0.0
    avg = 0
    for key, value in GToneFreq:
        sum = sum + value
        num = num +1
    if (sum > 0):
        avg = round (sum /num , 2) 
    return avg
def estimateGenomeSize(GToneFreq):
    freqThreshold = 2
    num = 0.0
    for key, value in GToneFreq:
        if (value > freqThreshold):
            num = num +1
    return num
def findCutoff(avgKmer , genomeSize):
    cutoff = minCutOff
    if (genomeSize < minGenomeSize):
        return minCutOff
    if (avgKmer> 1):
        cutoff = round( (avgKmer -1 )/ math.log(avgKmer),2)
    if (cutoff > maxCutoff):
        return maxCutoff
    if (cutoff < minCutOff):
        return minCutOff
    return cutoff
def extractKmers( readFileName):
    freq = {}
    GToneFreq = []
    with open(readFileName, "rU") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)
            i = 0
            if (len (seq) < kmerSize):
                continue
            while (i<=len (record.seq) - kmerSize):
                kmer= seq[i:i+kmerSize]
                kmerRC = RC(kmer)
                if (kmer > kmerRC ):
                    kmer = kmerRC
                try:
                    freq [kmer] = freq [kmer] + 1   
                except :
                    freq [kmer] = 1
                i = i +1
    for (key, value)  in freq.items():
        if (value > 1):
            t = (key, value)
            GToneFreq.append(t)
    return (GToneFreq)
def getBestKmerSize(cov ):
    if (cov < minCov):
        return minKmerSize
    coeff = math.floor( (cov - minCov)/6)
    newkmer = coeff * 5  + minKmerSize
    if (newkmer %2 == 0):
        newkmer = newkmer - 1
    if (newkmer > maxKmerSize):
        newkmer = maxKmerSize
    return newkmer
def main(argv=None):
    if argv == None:
        argv = sys.argv
    if len(argv) < 1:
        exit()
    readFile = str(sys.argv[1])
    #print ("read file :"+ readFile)
    outFile = str(sys.argv[2])
    #print ("directory :"+ outFile)
    trimReads(readFile,outFile)
    GToneFreq =  extractKmers(readFile)
    avg = computeAvg(GToneFreq)
    genomeSize=estimateGenomeSize(GToneFreq)    
    cutoff = findCutoff(avg, genomeSize)
    print (cutoff)
    kmer = getBestKmerSize(avg)
    print (kmer)  
    print (genomeSize)
main() 





 

