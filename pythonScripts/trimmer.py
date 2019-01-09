

#GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
#Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
#Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.


import sys
from Bio import SeqIO
from statistics import mean, median , pstdev
import os
import shutil
import math
kmerSize = 31
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
def findCutoff(avgKmer ):
    cutoff = round( (avgKmer -1 )/ math.log(avgKmer),2)
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
    GToneFreq =  extractKmers(outFile)
    avg = computeAvg(GToneFreq)    
    cutoff = findCutoff(avg)
    print (cutoff)

main() 





 
