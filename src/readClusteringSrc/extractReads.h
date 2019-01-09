
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

#include <vector>
#include <map>
#include <set>
#include <string>
#include <thread>
#include <mutex>
#include <algorithm>
#include <iomanip>
#include "../tkmer.h"
#include "../library.h"

using namespace std;
class fileManagment{
         mutex  fileMutex;                     // shared mutex between threads to make sure they write exclusively in the files 
         string tempDir;
         string normalFileName;
         string lowComplexFileName;
         ofstream normalFileStream;
         ofstream lowComplexFileStream;
public :
        fileManagment(LibraryContainer &libContRead , string tempDir_);
        void writeBackReads(vector<ReadRecord> &normalRcords,vector<ReadRecord> &lowComplex );
        ~fileManagment (){
                normalFileStream.close();
                lowComplexFileStream.close() ;
        };
};
class extractReads{
private:
        
        size_t kmerSize;                            // the given size of kmer    
        string normalFileName;                      // the normal file name, this thread wirts reads that dont' have any low comlex kmer here     
        string lowComplexFileName;                  // the low complex file name, this thread wirts reads that dont' have any low comlex kmer here         
        set <Kmer> kmerSet; 
        fileManagment &fileManager;                  // shared mutex between threads to make sure they write exclusively in the files 
        
        /**
         * check a pair of reads returns true if any of them has low complex kmer
         * @param forwardRead the first read from the pair
         * @param reverseRead the second read from the pair
         * @return true if they have low comlex kmer 
         **/
        bool  checkPairs(ReadRecord forwardRead , ReadRecord reverseRead);
        
        /**
         * check a single read to see if it has any lowComplex kemr 
         * @param read  the given read
         * @return true if it has low comlex kmer
         **/
        bool  checkRead (string read) ;
        void  writeNormalReads(vector<ReadRecord> &normal);
        void  wirteLowComplexReads(vector<ReadRecord> &lowComplex);
public :
        void checkReads (vector<ReadRecord> myReadBuf);
        extractReads (size_t kmer_ , set <Kmer > kmerSet_ , fileManagment &fileMa ) ;
        
};
class extractReadsHandler {
        size_t kmerSize;                        // the given size of kmer 
        string tempDir;                         // temporary directory to save the output files
        map <Kmer , size_t> kmerMap;            // shared dictonary, contains kmers and for the same kmers (i.e for the kmer and RC) returns a same value 
        string lowComplexFileName;              // the lowComplex file name , given to all threads to append their finding here
        string normalFileName;                  // the normal file name , given to all threads to append their finding here     
        std::mutex syncMutex;                   // shared mutex to make sure every time we read even number of reads for the pair reads       
        std::mutex fileMutex;                   // shared mutex between threads to make sure they write exclusively in the normal file
        set <Kmer> kmerSet; 
        /**
         * the files need to be in fasta format
         * load the given kmers from the file, kmers are assumed to be low complex like poly A
         * @param kmerFileName kmerFileName
         * @parma kmerList the output list, contains all the kmers
         **/
        void loadKmers(string kmerFileName, vector <string> &kmerList);
        /**
         * Extract the subkmers and RC of the given kmer list , for exame ATCATC is equal to TCATCA and CATCAT
         * @parma kmerList the input list, contains all the kmers
         **/
        void extractSubKmers(const vector <string> &kmerList);
        /**
         * extracts the reads that contains low complex kmers in parallel
         * @param libContSam the input library 
         * 
         * */
        void workerThread(LibraryContainer &libContSam , fileManagment &fileMa );
        /**
         * opens two output file, normal and lowComplex, later multiple threads append their reads into these files
         * the extension of these two files are based on the input format
         * @param libContRead contains the input read library
         * */
        void openFiles(LibraryContainer &libContRead);
        
        /**
         * reads the kmer file , and keeps the frequency of each kmer, we assume the frequency of each kmer is written in the id in file
         * @param kmerFileName kmerFileName
         * @param kmerFreq map keeps the frequency of each kmer 
         * */
        void loadKmersFreq(string kmerFileName ,  map <Kmer , unsigned int>& kmerFreq);
public :
        /**
        * the constructor
        * @param kmerFileName the name of kmer file, contains all the given kmer 
        * @param readFileName the name of read file, contains all the reads which needs to be examined if they contain any of those kmers
        * @param kmer the size of kmer
        * @param tempDir saves the output in this directory
        **/
        extractReadsHandler(string kmerFileName, string readFileName,size_t kmer , string tempDir);
        /**
         * manages the threads to extract reads that contain low complex kmer  from the file  
         * @param readFileName the name of read file, contains all the reads which needs to be examined if they contain any of those kmers
         **/
        void extractLowComplexReadsFromFile(string readFileName);
};
