

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
#include "extractReads.h"
#include <stdio.h>
#include "../util.h"
#include <pthread.h>
extractReadsHandler::extractReadsHandler(string kmerFileName, string readFileName, size_t kmer, string tmpDir):kmerSize(kmer), tempDir(tmpDir)
{
        Kmer::setWordSize(kmerSize);
        vector <string> kmerList;
        loadKmers(kmerFileName, kmerList);
        extractSubKmers(kmerList);
        extractLowComplexReadsFromFile(readFileName);
        
        
}
fileManagment::fileManagment(LibraryContainer& libContRead, string tempDir_):tempDir(tempDir_)
{
        const ReadLibrary &input = libContRead.getInput(0);
        std::ostringstream normalOss;
        normalOss << tempDir << "/normal." << input.getFileType() ;
        normalFileName =  normalOss.str();
        std::ostringstream lowComplexOss;
        lowComplexOss << tempDir << "/lowComplex." << input.getFileType() ;
        lowComplexFileName =  lowComplexOss.str();
        
        normalFileStream.open(normalFileName);
        lowComplexFileStream.open(lowComplexFileName) ;
}
void fileManagment::writeBackReads(vector<ReadRecord> &normalRcords,vector<ReadRecord> &lowComplexRecords ){
        std::lock_guard<std::mutex> lock(fileMutex);
        
        for (auto it:normalRcords){
                normalFileStream << it.preRead;
                normalFileStream << it.read;
                normalFileStream << it.postRead;
        }
        
        
        for (auto it:lowComplexRecords){
                lowComplexFileStream << it.preRead;
                lowComplexFileStream << it.read;
                lowComplexFileStream << it.postRead;
        }
        
        //fileMutex.unlock();
}

extractReads::extractReads (size_t kmer_  ,  set <Kmer > kmerSet_ , fileManagment &fileMa ) :kmerSize(kmer_), kmerSet(kmerSet_),  fileManager(fileMa)
{
        Kmer::setWordSize(kmerSize);

}

// we assume the file is in fasta format
void extractReadsHandler::loadKmers(string kmerFileName , vector <string> &kmerList)
{
        
        std::ifstream ifs(kmerFileName.c_str());
        size_t num = 0;
        while(!ifs.eof())
        {
                string line = "";
                string firstLine, secodnLine;
                std::getline(ifs, line);
                if (ifs.eof() || line == ""){
                        ifs.close();
                        break;
                }
                if (num %2 ==0){
                        firstLine = line ;
                        if (firstLine[0] != '>'){
                                cout << "Wrong format, the kmer file is not in fasta format " <<endl; 
                                exit(0);
                        }
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

                                cerr <<"kmer in the file is  smaller the given kmerSize  "<<kmerSize <<  "\n";
                                exit(0);
                        }
                        kmerList.push_back(line.substr(0, kmerSize));
                }
                num ++; 
        }
        ifs.close();
}



// we assume the file is in fasta format
void extractReadsHandler::loadKmersFreq(string kmerFileName ,  map <Kmer , unsigned int>& kmerFreq)
{
        
        std::ifstream ifs(kmerFileName.c_str());
        size_t num = 0;
        while(!ifs.eof())
        {
                string line = "";
                string firstLine, secodnLine;
                std::getline(ifs, line);
                unsigned int freq = 0;
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
                        string freqStr  = firstLine.substr(1, firstLine.length()-2);
                        cout << freqStr << endl;
                        
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
                        //kmerFreq.push_back(line);
                }
                num ++; 
        }
        ifs.close();
}


void extractReadsHandler:: extractSubKmers(const vector <string> &kmerList){
        
        for (auto kmerStr : kmerList){
                int i = 0 ;
                while (i <= kmerStr.length()){
                        char lastchar= kmerStr[kmerStr.length() -1];
                        kmerStr = lastchar + kmerStr.substr(0, kmerStr.length() -1);
                        Kmer kmer (kmerStr);
                        kmerSet.insert(kmer);
                        kmerSet.insert(kmer.getReverseComplement());
                        i = i+1;
                }
                /*for (auto it :kmerSet){
                        cout << it <<endl;
                } */   
                
        }
        
        /* for (std::string read :kmerList){
                // this will ensure that we can find patterns like this ATCGATCGATCG, TCGATCGATCGA,GATCGATCGATC
                for (KmerIt it(read); it.isValid(); it++){
                        Kmer kmer = it.getKmer();
                        kmerSet.insert(kmer);
                        kmerSet.insert(kmer.getReverseComplement());
                }
        }
        cout << "we are looking for these kmers in the data set. " <<endl;
        for (auto it :kmerSet){
                cout << it <<endl;
        }*/
}

void extractReads::checkReads (vector<ReadRecord> myReadBuf){
        assert(myReadBuf.size()%2 ==0 );
        vector<ReadRecord> lowComplex;
        vector<ReadRecord> normal;
        
        size_t i = 0;
        while (i < myReadBuf.size()){
                ReadRecord r1 = myReadBuf[i];
                ReadRecord r2 = myReadBuf[i+1];
   
                //write in low comlex file
                if (checkPairs(r1 , r2)){
                        lowComplex.push_back(r1);
                        lowComplex.push_back(r2);
                        
                }else{
                        normal.push_back(r1);
                        normal.push_back(r2); 
                }
                i = i +2;
        }
        fileManager.writeBackReads( normal, lowComplex);
        
}

bool extractReads::checkRead (string read ){
        if (read.size() < kmerSize)
                return false;
        for (KmerIt it(read); it.isValid(); it++){
                if (kmerSet.find(it.getKmer())!= kmerSet.end())
                        return true;
        }
        return false;
}
bool extractReads::checkPairs(ReadRecord r1 , ReadRecord r2){
        return (checkRead( r1.read) || checkRead( r2.read));
}
void extractReadsHandler::workerThread(LibraryContainer &libCon , fileManagment &fileMa ){
        extractReads extract(kmerSize, kmerSet, fileMa);
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
                extract.checkReads(myReadBuf);
        }
}
void extractReadsHandler::extractLowComplexReadsFromFile(string readFileName){
        size_t numOfThreads =  std::thread::hardware_concurrency();
        cout << "number of threads : " << numOfThreads <<endl;;
        LibraryContainer libContRead;
        libContRead.insert(ReadLibrary(readFileName,"",tempDir));
        size_t threadSize = 500000;
        libContRead.startIOThreads(threadSize, 10 * threadSize * numOfThreads, false);    
        // start worker threads
        vector<thread> workerThreads(numOfThreads);
        vector <ReadRecord> lastRecord;
        fileManagment fileMa(libContRead,tempDir );
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&extractReadsHandler::workerThread,this, ref(libContRead) , ref(fileMa));
        }
        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));
        libContRead.joinIOThreads();
        
}
