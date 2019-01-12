
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


#include<iostream>
#include "math.h"
#include"undGraph.h"
#include <assert.h>
#include <stack>
#include "../util.h"
#include <iomanip> 
#include <algorithm> 
#include <stdlib.h>

undGraph::undGraph( unsigned  int  numOfnodes)
{
        numNodes = numOfnodes + 1;
        adj = new vector< pair < unsigned  int , char >>[numNodes];
}
void undGraph::DFSItr( vector < set < unsigned  int  >> &components ,  unsigned  int  threshold)
{
        // Mark all the vertices as not visited
        bool *visited = new bool[numNodes];
        for( unsigned  int  v = 1; v < numNodes; v++)
                visited[v] = false;
        for ( unsigned  int  v = 1; v < numNodes; v++)
        {
                if (visited[v] == true)
                        continue;
                set < unsigned  int  > newComponent ;
                std::stack<unsigned int> stck;
                stck.push(v);
                while (!stck.empty()){
                        size_t source =  stck.top();
                        stck.pop();
                        visited[source] = true;
                        newComponent.insert(source);
                        vector <pair < unsigned  int , char>>::iterator it ;
                        for(it = adj[source].begin(); it != adj[source].end(); ++it){
                                if (it->second <= threshold) // weight = it->second;
                                        continue;
                                if (visited[it->first])      //target = it->first;
                                        continue;
                                stck.push (it->first);
                                visited[it->first] = true;
                        }
                }
                if (newComponent.size() > 0 )
                        components.push_back(newComponent);
        }
}

bool undGraph::addEdge( unsigned  int  source,  unsigned  int  target ,  float  weight)
{
        vector <pair< unsigned  int , char>> &currentEdges =  adj[source];
        bool found = false;
        for ( unsigned  int  i = 0; i < currentEdges.size() ; i++ ){
                pair< unsigned  int , char> &item = currentEdges [i];
                if (item.first == target ){
                        assert(item.second + weight < 255); // we use float 
                        found = true;
                        item.second = item.second + weight;
                        break;
                }
                
        }
        if (! found) { //this edge is new 
    
                adj[source].push_back(make_pair( target , weight));
                adj[target].push_back(make_pair( source , weight));
                return true;
                
        }else{     // edge is not new so increase the wight from the other direction as well
                vector <pair< unsigned  int , char>> &currentEdges =  adj[target];
                for ( unsigned  int  i = 0; i < currentEdges.size() ; i ++ ){
                        pair< unsigned  int , char> &item = currentEdges [i];
                        if (item.first == source ){
                                item.second = item.second + weight;
                                break;
                        }
                        
                }
                
        }
        return false;
}

void graphHandler::addEdgeThread(size_t threadID ,const ClusterReadAss &clusterReadMap, undGraph &g ){
        size_t edgesNum = 0;
        for (size_t i= clusterReadMap.getMinClusterID(); i<= clusterReadMap.getMaxClusterID(); i++){
                set <size_t > readSet = clusterReadMap.getReadsInCluster(i);
                std::set<size_t>::iterator firstIT, secondIT;
                for (firstIT = readSet.begin(); firstIT != readSet.end(); firstIT++)
                {
                        secondIT = firstIT;
                        secondIT ++;
                        while (secondIT != readSet.end()){
                                bool added = false;
                                if ( *firstIT >=  *secondIT ){
                                        if ((*firstIT % numThreads) == threadID ){
                                                if (g.addEdge(*firstIT, *secondIT, 1))
                                                        added = true;
                                        }
                                        
                                }else{
                                        if ((*secondIT % numThreads) == threadID ){
                                                if (g.addEdge(*secondIT, *firstIT, 1))
                                                        added = true;
                                        }
                                }
                                if (added){
                                        edgesNum ++;
                                        if (edgesNum % 1000 == 0 && threadID == 0 ){
                                                cout << "number of edges added to the graph: \t" << edgesNum <<"\r"  ;
                                                cout.flush();
                                        }
                                }
                                
                                secondIT++;
                        }
                }
        }
   
}
void graphHandler::buildGraph( vector <string> fileNames)
{
        Util::startChrono();
        cout << "number of threads : " <<numThreads <<endl;
        vector<thread> workerThreads(numThreads);
        vector<undGraph> graphArray(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++){
                undGraph g (G.getNumOfNodes());  
                graphArray [i] = g;
        }
        for (size_t i = 0 ; i < fileNames.size() ; i++){
                //map<size_t , set <size_t>> clusterReadMap;
                size_t numOfClusters  = getNumOfClusters(fileNames[i]);
                
                ClusterReadAss clusterReadDic(numOfClusters);
                readAdjFromFile(fileNames[i], clusterReadDic); 
                for (size_t i = 0; i < workerThreads.size(); i++){
                        workerThreads[i] = thread(&graphHandler::addEdgeThread, this,i ,ref (clusterReadDic), ref (graphArray[i]) );
                        
                }
                for (size_t i = 0; i < workerThreads.size(); i++){
                        workerThreads[i].join();        
                }
        }
        cout << "\nMerging edges from different threads. " <<endl;
        for (unsigned int source = 1; source <= G.getNumOfNodes(); source++){
                
                if (source % 1000 == 0  ){
                        
                        cout  <<"number of edges added to the graph: \t" <<std::fixed << std::setprecision(1)<<  ( (double)(source*100) /(double)G.getNumOfNodes() ) <<"\r"  ;
                        cout.flush();
                }
                
                for (size_t j = 0; j < workerThreads.size(); j++)
                        G.appendEdgeList(source, graphArray [j].getEdgeList(source));
                
        }
        cout << "\nBuilding graph finished successfully" <<endl;
}


size_t  graphHandler::getNumOfClusters( string  fileName){
        
        size_t maxClusterID = 1;
        std::ifstream ifs(fileName.c_str());
        map <size_t , set<size_t>> readCluster;
        while(!ifs.eof())
        {
                string line = "";
                
                std::getline(ifs, line);
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                vector <string > elements = Util::splitString(line , ' ');
                size_t clusterID = stoi(elements [1]);
                if (maxClusterID < clusterID){
                        maxClusterID = clusterID;
                }
        }
        ifs.close();  

        return maxClusterID;
}

void  graphHandler::readAdjFromFile( string  fileName, ClusterReadAss &clusterReadMap ){
        string firstLine = "" ;
        std::ifstream ifs(fileName.c_str());
        map <size_t , set<size_t>> readCluster;
        cout << "Reading from file " << fileName <<endl;
        while(!ifs.eof())
        {
                string line = "";
                
                std::getline(ifs, line);
                if (ifs.eof()){
                        ifs.close();
                        break;
                }
                vector <string > elements = Util::splitString(line , ' ');
                size_t  read = stoi(elements [0]);
                size_t clusterID = stoi(elements [1]);
                
                clusterReadMap.addreadToCluster(read,clusterID);
        }
        ifs.close();  

        cout << "Reading nodes and edges from  file finished !!  " <<endl;

}

void undGraph::print(string name )
{
        cout << "print Graph into disk file" <<endl;
        ofstream cytoGrpah( name);
        cytoGrpah << "source\ttarget\tweight" <<endl; 
        for ( unsigned  int  i = 0; i < numNodes; i++){
                vector <pair < unsigned  int , char>> edgeList = adj[i];
                
                for ( unsigned  int  j = 0; j < edgeList.size(); j++){
                        cytoGrpah << i<< "\t" << edgeList[j] .first <<"\t" <<int( edgeList[j] .second)  <<endl; ;
                }
        }
        
}
graphHandler::graphHandler(unsigned int numOfNodes_ , vector <string> fileNames  ) :numOfNodes(numOfNodes_)
{
        numThreads = thread::hardware_concurrency();
        G = undGraph (numOfNodes);
        buildGraph( fileNames);
        
}
int graphHandler::getOptimalClusterSize(double coverage, string staFileName)
{
        double avg_readLength = 0 ;
        size_t optimalClusterSize = 0;
        size_t kmerSize = 1;
        std::ifstream ifs(staFileName.c_str());
        std::set<string> ids;
        while(ifs.is_open() && !ifs.eof())
        {
                string line = "";
                std::getline(ifs, line);
                vector<string> items =Util::splitString(line,':');
                if(items[0]==  "AvgReadLen"){
                        avg_readLength = atof (items[1].c_str());
                }
                if(items[0]==  "kmerSize"){
                        kmerSize = atof (items[1].c_str());
                }
        }
        ifs.close(); 
        cout << "Avarage read length: " << avg_readLength <<endl;
        
         
        double errorEffect = pow(1-e , kmerSize );
        optimalClusterSize = round (errorEffect *(avg_readLength-kmerSize+1)*coverage/avg_readLength) ;
        cout << " The optimal cluster size is:\t" << optimalClusterSize <<endl;
        return optimalClusterSize;
}


void graphHandler:: extractClusters ( string outputClusteringFile, size_t iterNum , unsigned int optimalClusterSize  ){
        vector < set <unsigned int >> components ;
        components = G.getStableCore(iterNum, optimalClusterSize  , numOfNodes);
        std::map <unsigned int , unsigned int >   readCluster;
        unsigned int clusterID = 1;
        for (auto component:components){
                
                for (auto it :component){
                        readCluster[it] = clusterID;
                }
                clusterID ++;
        }
        cout << "Number of found clusters in the given interval size is : " << (clusterID-1) <<endl;
        ofstream output ( outputClusteringFile.c_str());
        for (unsigned int i = 1; i<= numOfNodes ; i++){
                if (readCluster.find(i) != readCluster.end()){
                        output << i << "\t" << readCluster.find(i)->second <<endl;;
                }
                else{
                        output << i << "\t" <<clusterID <<endl;;
                        clusterID ++;
                }
        }        
}


vector < set < unsigned  int  >> undGraph::getStableCore( size_t numOfClustering , unsigned int optimalClusterSize ,  unsigned  int  numNodes ){
        set< unsigned  int > visited ;
        vector < set < unsigned  int  >> finalComponents ;
        double maxMaxSize = 10 * (optimalClusterSize + 3*sqrt(optimalClusterSize)) ;//we assume at most 10 cluster are in one neighborhood
        double minMinSize =      (optimalClusterSize - 3*sqrt(optimalClusterSize)) > 10 ?(optimalClusterSize - 3*sqrt(optimalClusterSize)) : 10;  
        int minClusterSizeInitial = optimalClusterSize - 2*sqrt(optimalClusterSize);
        int maxClusterSizeInitial = optimalClusterSize + 2*sqrt(optimalClusterSize);
        double ratio = pow( (double)maxMaxSize / (double) maxClusterSizeInitial,(double)( 1/ (double)numOfClustering)) ;
        double  decrementRatio = (minClusterSizeInitial - minMinSize)/ (numOfClustering);
        unsigned  int  maxClusterSize =  maxClusterSizeInitial;
        int  minClusterSize = minClusterSizeInitial ;
        for (size_t i = 0; i < numOfClustering ; i++){
                vector < set < unsigned  int  >> components ;
            
                cout << "[" <<minClusterSize << ":" <<maxClusterSize << "]" <<endl;
                DFSItr (components, i );
                int  j = 0 ;
                for (auto component :components){
                        std::set< unsigned  int >::iterator it;
                        bool alreadyVisited = false;
                        for (it = component.begin(); it!= component.end(); it++)
                        {
                                if (visited.find(*it) != visited.end ()){
                                        alreadyVisited = true;
                                        break;
                                }
                        }
                        if (!alreadyVisited){
                                if ( component.size () <= maxClusterSize || i == numOfClustering-1){
                                        for (it = component.begin(); it!= component.end(); it++)
                                        {
                                                visited.insert(*it);
                                        }
                                }
                                if ( component.size () <= maxClusterSize && component.size() >= minClusterSize){
                                        finalComponents.push_back(component);
                                        j  = j +1;
                                }
                                if (i == numOfClustering-1 && component.size() > maxClusterSize){ // if this is the last round, report the biggest one
                                        
                                        finalComponents.push_back(component);
                                        j = j+1;
                                }
                        }
                }
                cout << "number of handeled Reads "<< visited.size()<< " / " << numNodes <<endl;
                cout << "number of connected component found at the level of " << i << " is, " << j <<endl;
                if (numNodes < visited.size() + minClusterSize ){ //early stop
                        break;
                }
                maxClusterSize =  maxClusterSizeInitial * pow (ratio, (i+1));
                minClusterSize = minClusterSizeInitial - decrementRatio * (i+1) ;
        }
        return finalComponents;
}


