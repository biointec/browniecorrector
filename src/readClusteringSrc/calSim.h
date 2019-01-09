


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

#include <mutex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>

#include <math.h>
#include <iomanip>
using namespace std;
const float  maxSimThreshold = 20;      
const float minSimThreshol = 10;
class ReadPair
{
public:
        string read1, read2;
        bool R1fwd, R1rev, R2fwd, R2rev;
        
public:
        /**
         * Default constructor
         */
        ReadPair() {};
        
        /**
         * Constructor
         */
        ReadPair(string read1, string read2) : read1(read1), read2(read2) {}
        
        void findPattern(const string& pattern) {
                string patternRC = revComp(pattern);
                
                R1fwd = (read1.find(pattern) != string::npos);
                R1rev = (read1.find(patternRC) != string::npos);
                R2fwd = (read2.find(pattern) != string::npos);
                R2rev = (read2.find(patternRC) != string::npos);
                
                if (!(R1fwd || R1rev || R2fwd || R2rev))
                        throw runtime_error("Some pair contains no pattern");
        }
        
        string revComp(const string& sequence)
        {
                string retval = sequence;
                reverse(retval.begin(), retval.end());
                
                for (size_t i = 0; i < retval.size(); i++) {
                        if (retval[i] == 'A')
                                retval[i] = 'T';
                        else if (retval[i] == 'T')
                                retval[i] = 'A';
                        else if (retval[i] == 'C')
                                retval[i] = 'G';
                        else if (retval[i] == 'G')
                                retval[i] = 'C';
                        else if (retval[i] == 'N')
                                retval[i] = 'N';
                        else
                                throw runtime_error("Trying to reverse complement non-ACGT sequences");
                }
                
                return retval;
        }
        
};

class Matrix
{
private:
        vector<int> matrix;
        int m;

public:
        /**
         * Constructor
         * @param m Number of rows
         * @param n Number of columns
         */
        Matrix(size_t m, size_t n) : matrix(m*n), m(m) {}

        /**
         * Operator () overloading
         * @param i Row index
         * @param j Column index
         * @return Element at position (i, j)
         */
        int operator() (int i, int j) const {
                return matrix[j*m + i];
        }

        /**
         * Operator () overloading
         * @param i Row index
         * @param j Column index
         * @return Reference to element at position (i, j)
         */
        int& operator() (int i, int j) {
                return matrix[j*m + i];
        }
};



class Alinger {
public : 
        Alinger(){};
        /**
         * Perform global alignment of two sequences and print the alignment to stdout
         * @param X sequence one
         * @param Y sequence two
         */
private :
        const int G = -7;
        const int M =  1;
        const int I = -4;
        int alignOverlap(const string& X, const string& Y);
public:
        int computeScore(const ReadPair& p1, const ReadPair& p2);

        
};

class calSim
{
        
private:
        size_t kmerSize;                                                //the kmer size which the graph is built with
        double initialCov;                                              // the initial coverage of the data set
        vector<size_t> nodeLenVec;                                      //vector contains lenght of nodes
        std::vector< vector<size_t> > readNodeSetRefined;               //for every read, keeps the associated valid nodes 
        std::vector< pair< vector<size_t>, vector<size_t>> > readNodeSetRefinedPair;  //for every read, keeps the associated valid nodes for each pair separately
        std::vector< vector < pair<size_t , size_t>> > readNodeSetEndPointsLen;   //for every read, for the last and first nodes keeps the aligned reads 
        vector<pair < pair<size_t, size_t> , double >> edgeList_all;    // the similarity score list, for each two reads there is one score
        size_t numThreads;                                              // number of threads 
        size_t numOfReads;                                              // number of reads (pairs) in the data set
        double coverage;
        double avg_readLength;
        size_t totalNumOfReads;
        size_t maxReadLen;
        vector<ReadPair> readPairs;
        Alinger aligner;
        std::map<size_t,size_t> kmerNodeCov;    
        
        /**
         * read the node file produced by brownieAligner and stores the lenght of nodes.
         * @nodeFile the node file name 
        
         **/
        void setNodeLen(string nodeFile);
        /**
         * according to the initial coverage , node lenght calculate the accepted range of coverage for each node
         * @param threshold  the accepted range of coverage is added to this value (for the lower bound is subtracted)
         * @param coverageVec for each lenght > minNodeLen keeps the accepted range of coverage 
         **/
        void setExpectedCov(float threshold , vector <pair < double , double> > & coverageVec);
        /**
         * based on the accepted coverage interval, for each read find the set of valid nodes that alined to
         * @param ncfFile the input alignment file, for each read keeps all the alined node
         * @param readNodeAssFile ouput file , to store the association of reads and ndoes (refined version of ncfFile)
         * 
         * */
        
        void  setReadNodeSet(string  ncfFile , std::vector< pair< set<size_t> , set<size_t>>>& readNodeSet );
        
        void  calEndpointLen(vector<string>  nodeList, int pos, vector <pair<size_t , size_t>> &endpoints , set <size_t> &nodeSet);
        void  setRefinedReadNodeSet(string ncfFile ,string readNodeAssFile ); 
        
        /**
         * finds the interction of two sets
         * @param v1    the first set
         * @param v2    the second set
         * */
        vector <size_t>  intersectionOfTwoSortedVe(vector <size_t> v1, vector <size_t> v2 );
        /**
         * finds the union  of two sets
         * @param v1    the first set
         * @param v2    the second set
         **/        
        vector <size_t > unionOfTwoSortedVe(vector <size_t> v1, vector <size_t> v2);
        
        /**
         * wirtes the association between nodes and reads in file
        
         * @param readNodeAssFile the file name which stores the read node associations
         * */
        void writeNodeReadAssInFile( string readNodeAssFile);
        /**
         * calculate the similarity score of the given read with all the rest read after itself
         * @param i the id of the current read
         * @param numOfReads the total number of reads
         **/
        vector<pair < pair<size_t, size_t> , double >> calSimForRead (size_t i , size_t numOfReads);
        vector<pair < pair<size_t, size_t> , double >> calAliSimForRead (size_t i , size_t numOfReads );
        double getSim(vector <size_t >& commonSet,vector<pair <size_t , size_t>>& firstEndPonts,vector<pair <size_t , size_t>>& secondEndPonts);

        /**
         * calculate the similarity score of those reads with other whose  read_id % numofthread == threadID
         * @param threadID the id of the current thread
         * @param mergeMutex mutex to integrate results from all threads
         */
        void calSimThread(size_t threadID ,  mutex & mergeMutex );
        
        void writeNetworkInFile(string networkFileName, size_t iterationNum=20);
        
        void uploadStatisticFromFile(string fileName);

        
public:
        
        calSim(string readFileName,string ncfFile_,  string nodeFile_, size_t kmerSize_, string readNodeAssFile_  , double initialCov_, string statFile);
        /**
         * calculate the similarity score of all to all reads
         * @param readNodeSetRefined for every read, keeps the associated valid ndoes
         * @param networkFileName output File name which stores the similarity scores for all binary combinations of reads
         * @param kmerSize the kmer size that graph built with
         * 
         * */
        void calSimInParallel( string networkFileName  );
        void readSequences(const string& filename, vector<ReadPair>& readPairs);
        
}; 
