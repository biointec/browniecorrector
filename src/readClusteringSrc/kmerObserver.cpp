
/************************************************************************************
*    Copyright (C) <2018> Mahdi Heydari <mahdi.heydari@ugent.be>                    *    
*                                                                                   *
*    This program is free software: you can redistribute it and/or modify           *
*    it under the terms of the GNU Affero General Public License as published       *
*    by the Free Software Foundation, either version 3 of the License, or           *
*    (at your option) any later version.                                            *    
*                                                                                   *    
*    This program is distributed in the hope that it will be useful,                *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *    
*    GNU Affero General Public License for more details.                            *
*                                                                                   *            
*    You should have received a copy of the GNU Affero General Public License       *
*    along with this program.  If not, see <https://www.gnu.org/licenses/>.         *
************************************************************************************/


#include <iostream>
#include <fstream>
#include "kmerObserver.h"
#include <stdio.h>
#include <pthread.h>
extractKmersHandler::extractKmersHandler(string kmerFileName, string readFileName, size_t kmer, string tmpDir):kmerSize(kmer), tempDir(tmpDir)
{
        Kmer::setWordSize(kmerSize);

        loadKmersFreq( kmerFileName);
        extractKmerFreqFromFile(readFileName);
        
        
}

extractKmers::extractKmers (size_t kmer_  , map <Kmer , unsigned int> kmerFreq_in  ) :kmerSize(kmer_) , kmerFreq(kmerFreq_in)
{
        Kmer::setWordSize(kmerSize);
}



// we assume the file is in fasta format
void extractKmersHandler::loadKmersFreq(string kmerFileName)
{
        
        std::ifstream ifs(kmerFileName.c_str());
        size_t num = 0;
        unsigned int freq = 0;
        while(!ifs.eof())
        {
                string line = "";
                string firstLine, secodnLine;
                std::getline(ifs, line);
                
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                if (num %2 ==0){
                        firstLine = line ;
                        if (firstLine[0] != '>'){
                                cout << "Wrong format, the kmer file is not in fasta format " <<endl; 
                                exit(0);
                        }
                        string freqStr  = firstLine.substr(1, firstLine.length()-1);
  
                        freq = stoi(freqStr);
                        
                }
                else{
                        secodnLine = line;
                        for (size_t i = 0; i< line.length(); i++){
                                if (toupper(line[i])!='A' && toupper(line[i])!='T'
                                        && toupper(line[i])!='C' && toupper(line[i])!='G' ){
                                        cout << "Error occurred, the sequence should contain only A, T, C, G "  <<endl;
                                exit(0);
                                        }
                        }
                        if (line.size() < kmerSize){
                                cout << "Error occurred, the sequence length should be at least as large as kmer size : " << kmerSize  <<endl;
                                exit(0);
                        }
                        string  kmer = line;
                        string  kmerRC = kmer;
                        Nucleotide::revCompl(kmerRC);
                        if (kmer < kmerRC)
                                kmerFreq [kmer] = freq;
                        else
                                kmerFreq [kmerRC] = freq;        
                        freq = 0;
                }
                num ++; 
        }
        ifs.close();
}



void extractKmers::checkReads (vector<ReadRecord> &myReadBuf){

        
        size_t i = 0;
        while (i < myReadBuf.size()){
                checkRead(myReadBuf[i].read);
                i ++;
        }
        
}

bool extractKmers::checkRead (string & read ){
        if (read.size() < kmerSize)
                return false;
        string freqStr = "";
        for (KmerIt it(read); it.isValid(); it++){
                Kmer kmer = it.getKmer();
                Kmer kmerRC = kmer.getReverseComplement();
                unsigned int freq = 1;
                if (kmer < kmerRC){
                        if (kmerFreq.find (kmer) != kmerFreq.end())
                                freq = kmerFreq[kmer];// [kmer] ;
                        
                }
                else{
                        if (kmerFreq.find(kmerRC) != kmerFreq.end())
                                freq = kmerFreq[kmerRC]; 
                        
                }
                freqStr = freqStr + " " + to_string( freq);
        }
        //cout << freqStr <<endl;
        read = freqStr ;
        return true;
}

void extractKmersHandler::workerThread(LibraryContainer &libCon ){
        extractKmers extract(kmerSize,  kmerFreq);
        bool result = true;
        while (result) {
                size_t blockID, recordID;
                vector<ReadRecord> myReadBuf;
                
                result = libCon.getRecordChunk(myReadBuf, blockID, recordID);
                if (myReadBuf.size() % 2) {
                        std::cout << "buf size " << myReadBuf.size()<<std::endl;
                }
                
                
                if (!result )
                        break;
                else{
                        extract.checkReads(myReadBuf);
                        libCon.commitRecordChunk(myReadBuf, blockID, recordID);
                }
                
        }
}
void extractKmersHandler::extractKmerFreqFromFile(string readFileName){
        size_t numOfThreads =  1 ; //std::thread::hardware_concurrency();
        cout << "number of threads : " << numOfThreads <<endl;;
        LibraryContainer libContRead;
        libContRead.insert(ReadLibrary(readFileName,"",tempDir));
        size_t threadSize = 500000;
        libContRead.startIOThreads(threadSize, 10 * threadSize * numOfThreads, true);    
        // start worker threads
        vector<thread> workerThreads(numOfThreads);
        vector <ReadRecord> lastRecord;

        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&extractKmersHandler::workerThread,this, ref(libContRead) );
        }
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        libContRead.joinIOThreads();
        
}
