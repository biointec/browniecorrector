
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

class extractKmers{
private:
        map <Kmer , unsigned int> kmerFreq ;
        size_t kmerSize;                            // the given size of kmer    

        
        /**
         * check a single read to see if it has any lowComplex kemr 
         * @param read  the given read
         * @return true if it has low comlex kmer
         **/
        bool  checkRead (string & read) ;

public :
        void checkReads (vector<ReadRecord> &myReadBuf);
        extractKmers (size_t kmer_ , map <Kmer , unsigned int> kmerFreq_in ) ;
        
};
class extractKmersHandler {
        size_t kmerSize;                        // the given size of kmer 
        string tempDir;                         // temporary directory to save the output files
        map <Kmer , unsigned int> kmerFreq ;

        /**
         * extracts the reads that contains low complex kmers in parallel
         * @param libContSam the input library 
         * 
         * */
        void workerThread(LibraryContainer &libContSam  );

        
        /**
         * reads the kmer file , and keeps the frequency of each kmer, we assume the frequency of each kmer is written in the id in file
         * @param kmerFileName kmerFileName
         * @param kmerFreq map keeps the frequency of each kmer 
         * */
        void loadKmersFreq(string kmerFileName );
public :
        /**
        * the constructor
        * @param kmerFileName the name of kmer file, contains all the given kmer 
        * @param readFileName the name of read file, contains all the reads which needs to be examined if they contain any of those kmers
        * @param kmer the size of kmer
        * @param tempDir saves the output in this directory
        **/
        extractKmersHandler(string kmerFileName, string readFileName,size_t kmer , string tempDir);
        /**
         * manages the threads to extract reads that contain low complex kmer  from the file  
         * @param readFileName the name of read file, contains all the reads which needs to be examined if they contain any of those kmers
         **/
        void extractKmerFreqFromFile(string readFileName);
};
