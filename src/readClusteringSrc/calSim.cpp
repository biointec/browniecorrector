
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
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <iomanip> 
#include "calSim.h"
#include <mutex>
#include <thread>
#include <pthread.h>
#include "../util.h"
#include <assert.h>

using namespace std;


inline double round( double val )
{
        if( val < 0 ) return ceil(val - 0.5);
        return floor(val + 0.5);
}
inline string findBetween(string str, string firstStr, string lastStr){
        
        unsigned first = str.find(firstStr);
        unsigned last = str.find(lastStr);
        string strBet = str.substr (first+ firstStr.length(),last-first -lastStr.length());
        return strBet;
}
calSim::calSim( string readFilename, string ncfFile_,string nodeFile_, size_t kmerSize_, string readNodeAssFile_  , double initialCov_,string staFile):kmerSize(kmerSize_), initialCov(initialCov_), avg_readLength(100), maxReadLen(100){
        numThreads = std::thread::hardware_concurrency();
       
        setNodeLen(nodeFile_);
        setRefinedReadNodeSet( ncfFile_, readNodeAssFile_);
        numOfReads = readNodeSetRefined.size() ;
        uploadStatisticFromFile(staFile);
        readSequences(readFilename, readPairs);
}
void calSim::setNodeLen(string nodeFile){
        
        std::ifstream ifs(nodeFile.c_str());
        size_t lineNum = 0;
        while(!ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                if (lineNum % 2 == 0){
                        vector<string> firstList =Util::splitString(line, '\t');
                        double cov = stof (firstList[3])/(stof(firstList[2]) -kmerSize +1);
                        kmerNodeCov [stoi(firstList[1])] = cov;
                        nodeLenVec.push_back(stof(firstList[2]));
                }
                
                
                lineNum ++; 
        }
        ifs.close();
}
void calSim:: setReadNodeSet(string  ncfFile , std::vector< pair< set<size_t> , set<size_t>>>& readNodeSet ){
        
        string firstLine = "" ;
        string secodnLine = "" ;
        std::ifstream ifs(ncfFile.c_str());
        size_t num = 0;
        while(!ifs.eof())
        {
                string line = "";
                
                std::getline(ifs, line);
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                if (num %2 ==0)
                        firstLine = line ;
                else{
                        vector <pair<size_t , size_t>> endpoints; 
                        secodnLine = line;
                        string firstStr   = findBetween(firstLine, "C(", ")C");
                        string secondStr  = findBetween(secodnLine, "C(", ")C");
                        vector<string> firstList =Util::splitString(firstStr, ' ');
                        vector<string> secondList = Util::splitString(secondStr, ' ');
                        string firstPos   = findBetween(firstLine, "P(", ")P");
                        string secondPos  = findBetween(secodnLine, "P(", ")P");
                        vector<string> firstPosList  = Util::splitString(firstPos, ',');
                        vector<string> secondPosList = Util::splitString(secondPos, ',');
                        int pos1 = stoi( firstPosList[0]);
                        int pos2 = stoi( secondPosList[0]);
                        set <size_t> firstNodeSet ;
                        set <size_t> secondNodeSet ;
                        calEndpointLen(firstList, pos1,endpoints , firstNodeSet);
                        calEndpointLen(secondList, pos2, endpoints , secondNodeSet);
                        readNodeSet.push_back( make_pair( firstNodeSet, secondNodeSet));
                        readNodeSetEndPointsLen.push_back(endpoints);
                }
                num ++; 
        }
        ifs.close();
        
}

void calSim::calEndpointLen(vector<string>  nodeList, int pos, vector <pair<size_t , size_t>> &endpoints ,
          set <size_t> &nodeSet){
        size_t nodeNum = 0;
        int sumLen = 0;
        for (auto it : nodeList){
                int nodeId = abs( stoi(it));
                if (nodeId ==0){
                        nodeNum ++;
                        continue;
                }
                nodeSet.insert(nodeId);
                // if there is only one node in the chain
                if (nodeList.size() ==1){
                        int maxOverlap = avg_readLength-kmerSize+1  ;
                        assert(maxOverlap>0);
                        endpoints.push_back(make_pair(abs(nodeId),maxOverlap));
                        break;
                }
                // if there are more than one node in the chain
                int nodeLength = nodeLenVec[nodeId-1];
                if (nodeNum ==0){
                        int maxOverlap = nodeLength - pos- kmerSize+1 ;
                        assert(maxOverlap>0);
                        endpoints.push_back(make_pair(abs(nodeId),maxOverlap));
                        sumLen = nodeLength - pos;
                }
                if (nodeNum != 0 && nodeNum != nodeList.size()-1){
                        sumLen = sumLen + nodeLength -kmerSize+1;
                }
                if (nodeNum == nodeList.size()-1){
                        int maxOverlap = avg_readLength-sumLen;
                        assert(maxOverlap>0);
                        endpoints.push_back(make_pair(abs(nodeId),maxOverlap));
                }
                nodeNum = nodeNum +1;
        }
}
void calSim::writeNodeReadAssInFile(string readNodeAssFile){
        ofstream tableFile(readNodeAssFile.c_str());
        size_t i = 0;
        
        while (i < readNodeSetRefined.size()){
                vector <size_t > nodeSet = readNodeSetRefined[i];
                string line = to_string(i +1) + "\t";
                
                for (size_t  it : nodeSet){
                        string nodeStr =  to_string(it);
                        line = line +  nodeStr + " ";
                }
                tableFile <<line<<endl;
                i  = i +1;
                
        }
        tableFile.close();
}
void calSim:: setRefinedReadNodeSet(string ncfFile,string readNodeAssFile ){
        
         
        double  errorEffect = pow(1-errorRate , kmerSize );
        coverage = round (errorEffect *(avg_readLength-kmerSize+1)*initialCov/avg_readLength) ;
        double c_k = round (errorEffect *(avg_readLength-kmerSize+1)*coverage/avg_readLength);
        double minCov = c_k/2 - sqrt(c_k/2)*3 < 2 ? 2: c_k/2- sqrt(c_k/2)*3; 
        double maxCov = c_k/2 + sqrt(c_k/2)*3 > c_k ? c_k: c_k/2 + sqrt(c_k/2)*3;;
        std::vector< pair< set<size_t> , set<size_t>>> readNodeSet ;
        setReadNodeSet( ncfFile,  readNodeSet);
        while (maxCov <= c_k){
              
                readNodeSetRefined.clear();
                
                cout << "[min, max] : [ " <<minCov <<", " <<maxCov <<" ]\n";
                size_t notSelectedReads = 0;
                size_t notSelectedPairs = 0;
                size_t readIndex = 1;
                set <size_t> selectedNodesGlobal;
                for (auto it : readNodeSet){
                        pair <set <size_t> , set <size_t> >  nodeSet = (readNodeSet[readIndex -1]);
                        set <size_t>  selectedNodesFirPair;
                        set <size_t>  selectedNodesSecPair;
                        set <size_t>  selectedNodesLocal;
                        bool selected = false;
                        bool selectedPair = false;
                        for (auto node :nodeSet.first){
                                if (kmerNodeCov[node] > minCov && kmerNodeCov[node] < maxCov ){
                                        selectedNodesLocal.insert(node);
                                        selectedNodesGlobal.insert(node);
                                        selectedNodesFirPair.insert(node);
                                        selected =true;
                                }
                        }
                        if (selected == false){
                                notSelectedReads = notSelectedReads +1;
                        }
                        else{
                                selectedPair = true;
                        }
                        selected = false;
                        for (auto node :nodeSet.second){
                                if (kmerNodeCov[node] > minCov && kmerNodeCov[node] < maxCov ){
                                        selectedNodesLocal.insert(node);
                                        selectedNodesGlobal.insert(node);
                                        selectedNodesSecPair.insert(node);
                                        selected =true;
                                }
                        }
                        if (selected == false){
                                notSelectedReads = notSelectedReads +1;
                        }else{
                                selectedPair = true;
                                
                        }
                        if (!selectedPair )
                                notSelectedPairs = notSelectedPairs +1;
                        std::vector<size_t> selectedNode(selectedNodesLocal.size());
                        std::copy(selectedNodesLocal.begin(), selectedNodesLocal.end(), selectedNode.begin());
                        std::sort(selectedNode.begin(), selectedNode.end());
                        readNodeSetRefined.push_back (selectedNode   );
                         
                        std::vector<size_t> selectedNodeFirst(selectedNodesFirPair.size());
                        std::copy(selectedNodesFirPair.begin(), selectedNodesFirPair.end(), selectedNodeFirst.begin());
                        std::sort (selectedNodeFirst.begin(),selectedNodeFirst.end());
                        
                        std::vector<size_t> selectedNodeSecond(selectedNodesSecPair.size());
                        std::copy(selectedNodesSecPair.begin(), selectedNodesSecPair.end(), selectedNodeSecond.begin());
                        std::sort (selectedNodeSecond.begin(), selectedNodeSecond.end ());
                        readNodeSetRefinedPair.push_back (make_pair(selectedNodeFirst,selectedNodeSecond)  );
                        
                        
                       
                        readIndex ++;                
                }
                double ratioSingle =  double(notSelectedReads*100)/ (double)(2*(readIndex-1));
                double ratioPiar   =  double(notSelectedPairs*100)/ double((readIndex-1));
                maxCov ++;
                
                cout << "number of reads without any selected node : " <<  notSelectedReads << " out of " <<( 2*(readIndex-1)) << " ( "<< ratioSingle << " %)"<<  endl;
                cout << "number of pairs without any selected node : " <<  notSelectedPairs << " out of "  <<  (readIndex-1) << " ( " << ratioPiar << " %)"<< endl;
                cout << "number of selected nodes :"          << selectedNodesGlobal.size() <<endl ;
                if (ratioSingle < 2 || ratioPiar < 1   )
                        break;
                else{
                        cout << "There are too many reads uncovered, we need to increase the coverage interval !!" <<endl;
                }
        }
        writeNodeReadAssInFile(readNodeAssFile);
        cout <<endl;
}


inline vector <size_t > calSim::unionOfTwoSortedVe(vector <size_t> first, vector <size_t> second){
        
        std::vector<size_t> result;
        
        std::set_union(first.begin(), first.end(),
                       second.begin(), second.end(),                  
                       std::back_inserter(result));
        return result;
}

inline vector <size_t> calSim::intersectionOfTwoSortedVe(vector <size_t> v1, vector <size_t> v2 ){
        
        
        std::vector<size_t> v_intersection;
        
        std::set_intersection(v1.begin(), v1.end(),
                              v2.begin(), v2.end(),
                              std::back_inserter(v_intersection));
        return v_intersection;
}


void calSim::calSimInParallel(string networkFileName)

{      
        cout << "Finding the similarity score of all to all reads." <<endl;
        cout << "This might take a few minutes be patient ..." <<endl;
        Util::startChrono();
        std::mutex mergeMutex;
        cout << "number of threads : " <<numThreads <<endl;
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i] = thread(&calSim::calSimThread, this, i ,ref(mergeMutex) );
                
        }
        for (size_t i = 0; i < workerThreads.size(); i++){
                workerThreads[i].join();        
        }
        sort(edgeList_all.begin(), edgeList_all.end());

        writeNetworkInFile (networkFileName);
        cout << "Done similarity calculation ("
        << Util::stopChronoStr() << ")" << endl;
}

void calSim::writeNetworkInFile (string networkFileName, size_t iterationNum ){
        size_t pos = networkFileName.find_last_of("/\\");
        string base_filename = networkFileName.substr(pos+ 1);
        string path = networkFileName.substr(0,pos+1);
        pos = base_filename.find_last_of('.');
        string file_without_extension = base_filename.substr(0, pos);
        string extension = base_filename.substr(pos, base_filename.length()-1);
        

        for (size_t i = 0; i< iterationNum ; i++){
                networkFileName = path + file_without_extension + "_" + to_string(i) + extension;
                ofstream networkFile(networkFileName);
                for( auto edge:edgeList_all)
                {

                        if ( edge.second <=  (minSimThreshol + ((i*(maxSimThreshold-minSimThreshol))/iterationNum)) ){
                                continue;
                                
                        }
                        
                        networkFile << edge.first.first <<  "\t"  <<  
                        edge.first.second <<  "\t" <<
                        std::setprecision(2)<< std::fixed <<   edge.second << endl;
                } 
                networkFile.close(); 
        }
}
void calSim::calSimThread(size_t threadID, std::mutex &mergeMutex ){
        
        vector<pair < pair<size_t, size_t> , double >> edgeList_t; 
        for  (size_t i = 1; i <= numOfReads; i++){
                if (i% numThreads !=threadID )
                        continue;
                vector<pair < pair<size_t, size_t> , double >> edgeList_i = calAliSimForRead( i, numOfReads);
                edgeList_t.insert(edgeList_t.end(), edgeList_i.begin(),edgeList_i.end());
                
                if (i% 100 == 0)
                {
                        double perc = 100.0 * (double)i / (double)numOfReads;
                        cout << std::fixed << std::setprecision(1) << "\tProcessing file (" << perc << "%)\r";
                        cout.flush();
                }
        }
        
        std::unique_lock<std::mutex> lock(mergeMutex);
        edgeList_all.insert(edgeList_all.end(), edgeList_t.begin(),edgeList_t.end());
        lock.unlock();
        
}

vector<pair < pair<size_t, size_t> , double >> calSim::calAliSimForRead (size_t i , size_t numOfReads ){
        
        size_t j = i+1;
        vector<pair < pair<size_t, size_t> , double >> edgeList_i ;
        vector <size_t > firstSet = readNodeSetRefined[i-1];
        if ( firstSet.size() == 0 )
                return edgeList_i;
        ReadPair R1 = readPairs[i-1];
        while (j <= numOfReads){
                vector <size_t > secondSet = readNodeSetRefined[j-1];
                if (secondSet.size() == 0 ){
                        j = j + 1;
                        continue;
                }
                vector <size_t > intersection  = intersectionOfTwoSortedVe(firstSet, secondSet);
                if (intersection.size() == 0 ){
                        j = j + 1 ;
                        continue;
                }
                double sim =0;
                double maxSim = avg_readLength *2;
                ReadPair R2 = readPairs[j-1];
                sim = aligner.computeScore(R1, R2);
                //double ratio = maxSim/ (2*avg_readLength - 2*kmerSize+2);
                sim =  ( double (sim) / double(maxSim)) ;
                sim = sim *100;
                if (sim > 0)
                        edgeList_i.push_back(make_pair( make_pair(i, j), sim));
                j = j+1;
        }
        return edgeList_i;
}
void calSim::readSequences(const string& filename, vector<ReadPair>& readPairs)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        string qname1, read1, qname2, read2, dummy;
        while (ifs) {
                getline(ifs, qname1);     // qname 1
                if (qname1.empty())
                        continue;
                if (qname1.front() != '@')
                        throw runtime_error("File: " + filename + " is not an interleaved FASTQ file");
                getline(ifs, read1);     // read 1
                transform(read1.begin(), read1.end(), read1.begin(), ::toupper);
                getline(ifs, dummy);     // + sign
                getline(ifs, dummy);     // quality scores 1

                getline(ifs, qname2);     // qname 1
                if (qname2.front() != '@')
                        throw runtime_error("File: " + filename + " is not an interleaved FASTQ file");
                getline(ifs, read2);     // read 2
                transform(read2.begin(), read2.end(), read2.begin(), ::toupper);
                getline(ifs, dummy);     // + sign
                getline(ifs, dummy);     // quality scores 2
                if (qname1 != qname2)
                        throw runtime_error("File: " + filename + " is not an interleaved FASTQ file");
                readPairs.push_back(ReadPair(read1, read2));
        }
}



void calSim::uploadStatisticFromFile(string fileName){
        std::ifstream ifs(fileName.c_str());
        std::set<string> ids;
        while(ifs.is_open() && !ifs.eof())
        {
                string line = "";
                
                std::getline(ifs, line);
                vector<string> items = Util::splitString(line,':');
                if(items[0]==  "NumOfReads"){
                        totalNumOfReads = atoi (items[1].c_str())/2;
                }
                if(items[0]==  "AvgReadLen"){
                        avg_readLength = atoi (items[1].c_str());
                }
                if(items[0]==  "MaxReadLen"){
                        maxReadLen = atoi (items[1].c_str());
                }
                
        }
        cout << "Number of reads: " << totalNumOfReads <<endl;
        cout << "Avarage read length: " << avg_readLength <<endl;
        cout << "Max read length: " <<maxReadLen <<endl;
        ifs.close();
}

int Alinger::computeScore(const ReadPair& p1, const ReadPair& p2)
{
        // forward forward alignment
        int FF1 = alignOverlap(p1.read1, p2.read1) + alignOverlap(p1.read2, p2.read2);
        int FF2 = alignOverlap(p1.read2, p2.read1);
        int FF3 = alignOverlap(p1.read1, p2.read2);
        
        int maxFF = max<int>(FF1, max<int>(FF2, FF3));
        
        ReadPair p1RC;
        p1RC.read1 = p1RC.revComp(p1.read2);
        p1RC.read2 = p1RC.revComp(p1.read1);
        
        // forward forward alignment
        int RF1 = alignOverlap(p1RC.read1, p2.read1) + alignOverlap(p1RC.read2, p2.read2);
        int RF2 = alignOverlap(p1RC.read2, p2.read1);
        int RF3 = alignOverlap(p1RC.read1, p2.read2);
        
        int maxRF = max<int>(RF1, max<int>(RF2, RF3));
        
        return max<int>(maxFF, maxRF);
}
int Alinger::alignOverlap(const string& X, const string& Y)
{
  

        int m = (int)X.length();
        int n = (int)Y.length();

        // initialize an (m+1) x (n+1) matrix S
        Matrix S(m+1, n+1);

        // initialize first column
        for (int i = 0; i <= m; i++)
                S(i, 0) = 0;

        // initialize first row
        for (int j = 1; j <= n; j++)
                S(0, j) = 0;

        // fill in the rest of the matrix
        int maxVal = 0, maxi = m, maxj = 0;
        for (int i = 1; i <= m; i++) {
                for (int j = 1; j <= n; j++) {
                        int diag = S(i-1, j-1) + (X[i-1] == Y[j-1] ? M : I);
                        int gapX = S(i, j-1) + G;
                        int gapY = S(i-1, j) + G;
                        S(i, j) = max(max(diag, gapX), gapY);
                }

                if (S(i, n) > maxVal) {
                        maxVal = S(i, n);
                        maxi = i;
                        maxj = n;
                }
        }

        for (int j = 0; j <= n; j++) {
                if (S(m, j) > maxVal) {
                        maxVal = S(m, j);
                        maxi = m;
                        maxj = j;
                }
        }

        // create an alignment
        /*string alX, alY, mid;

        int i = maxi;
        int j = maxj;

        while (i > 0 && j > 0) {
                if ((i > 0) && (S(i, j) == S(i-1, j) + G)) {
                        alX.push_back(X[i-1]);
                        alY.push_back('-');
                        mid.push_back(' ');
                        i--;
                } else if ((j > 0) && (S(i, j) == S(i, j-1) + G)) {
                        alX.push_back('-');
                        alY.push_back(Y[j-1]);
                        mid.push_back(' ');
                        j--;
                } else {
                        alX.push_back(X[i-1]);
                        alY.push_back(Y[j-1]);
                        char c = (X[i-1] == Y[j-1]) ? '|' : '*';
                        mid.push_back(c);
                        i--;
                        j--;
                }
        }

        for (int c = 0; c < i; c++) {
                alY.push_back(' ');
                mid.push_back(' ');
        }

        for (int c = 0; c < j; c++) {
                alX.push_back(' ');
                mid.push_back(' ');
        }

        reverse(alX.begin(), alX.end());
        reverse(alY.begin(), alY.end());
        reverse(mid.begin(), mid.end());

        alX = X.substr(0, i) + alX + X.substr(maxi);
        alY = Y.substr(0, j) + alY + Y.substr(maxj);

        cout << "Overlap alignment X[" << i << "-" << maxi-1 << "], Y["
             << j << "-" << maxj-1 << "]\n";
        for (size_t i = 0; i < alX.size(); i += 250) {
                cout << alX.substr(i, 250) << "\n"
                     << mid.substr(i, 250) << "\n"
                     << alY.substr(i, 250) << "\n\n";
        }

        cout << "Alignment score: " << S(maxi, maxj) << endl;*/

        return S(maxi, maxj);
}
