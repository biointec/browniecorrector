
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

#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <thread>
#include <pthread.h>
// undGraph class represents a undirected graph
// using adjacency list representation

using namespace std;

class ClusterReadAss {
private :
        size_t maxNumClusters;
        set<size_t> * clusterRead;
        size_t maxCurClusterID;
        size_t minCurClusterID;
        size_t numOfAddedReads;
        set<size_t> clusterSet;
public:
        // We start clusterID from 1
        ClusterReadAss(size_t _maxNumClusters):maxCurClusterID(1),minCurClusterID(1),numOfAddedReads(0){
                maxNumClusters = _maxNumClusters + 1; 
                clusterRead = new set<size_t> [maxNumClusters];
        }
        ~ClusterReadAss(){
                delete [] clusterRead;
                clusterRead = NULL;
        };
        /**
         * add read ID to the clusterSet
         * @param readID the read ID 
         * @param clusterID the cluster ID
         */
        void addreadToCluster(size_t readID, size_t clusterID){
         
                //cout << readID << "\t"<<clusterID <<endl;
                clusterRead[clusterID].insert(readID);
                numOfAddedReads = numOfAddedReads +1;
                
                if (clusterID > maxCurClusterID){
                        maxCurClusterID = clusterID;
                }
                if (numOfAddedReads ==1){ //only do it for the first read
                        minCurClusterID = clusterID;
                }
                else{
                        if (clusterID <minCurClusterID )
                                minCurClusterID = clusterID;
                }
        }
        set<size_t>  getReadsInCluster(size_t clusterID)const {
                set<size_t>  result;
                if (clusterID >= minCurClusterID && clusterID <=maxCurClusterID){
                        result =  clusterRead[clusterID];
                }
                else{
                        std::cerr<<"Error occured, the index is out of range";
                }
                return result;
        }
        size_t getMinClusterID ()const{
                return minCurClusterID;
        }
        size_t getMaxClusterID ()const{
                return maxCurClusterID;
        }
        const set<size_t> * begin() const{
                return &clusterRead[minCurClusterID];
        }
        const set<size_t> * end() const{
                return &clusterRead[maxCurClusterID+1];
        }
        set<size_t> * operator []  (size_t clusterID)const{
                if (clusterID >= minCurClusterID && clusterID <=maxCurClusterID){
                        return  &clusterRead[clusterID];
                }
                else{
                        std::cerr<<"Error occured, the index is out of range";
                        return NULL;
                }
        }
        
};


class undGraph
{
         
        unsigned  int  numNodes;    // No. of nodes (paired reads) in the graph
        vector<pair < unsigned  int  , char >> *adj; //the adjacency list, contains the edges with weights, weight is max (255)
        set <unsigned int> subnetworkSet;
        
        /**
         * recursively traverse the graph to find the connected component from the source node
         * @param source the search starts from here
         * @param visited a pointer to an array contains all visited nodes , hence should be escaped
         * @param component keeps the nodes which are in this component 
         * @param threshold the threshold for escaping edges , if the weight is less than this value
         **/
        void DFS( unsigned  int  source, bool visited[],  set < unsigned  int > & component,  unsigned  int  threshold);



public:
        /**
         * constructor 
         * @param numOfNodes number of nodes in the grpah
         *
         * */
        undGraph( unsigned  int  numOfNodes);   // Constructor
        /**
         * default constructor 
         * @param numOfNodes number of nodes in the grpah
         *
         * */
        
        undGraph() {}
        /**
         * gives the min weight vlaue which can split the graph into two parts
         * 
         * */
               /**
         * add edge to graph by reading from file
         * @parma networkFile the input file contains all the edges with weight
         * 
         * */
         void addEdgesFromFile(string networkFile);
        /**
         * iteratively traverse the graph to find the connected component from the source node (non recursive function)
         * @param components keeps all the components
         * @param threshold minimum weight between nodes in one component
         **/        
        void DFSItr(  vector < set < unsigned  int  >> &components ,  unsigned  int  threshold);
  
        /**
         *  add a directed edge to the graph. the order of  source and target in the arguments are important.
         * @param source the source nodes
         * @param target the target nodes
         * @weight the weight of edge between two nodes
         * @return true if it adds a new edge in graph
         * */
        bool addDirEdge( unsigned  int  source,  unsigned  int  target ,  unsigned  int  weight);
        /**
         * add edge to the graph, it is not directed,  if it exist increase the weight for both direction
         * @param source the source nodes
         * @param target the target nodes
         * @weight the weight of edge between two nodes
         * @return true if it adds a new edge in graph
         * */
        bool addEdge( unsigned  int  source,  unsigned  int  target,  float  weight);
        /**
         * finds all connected component , all the edges between nodes have weight >= threshold
         * @param components keeps all the components
         * @param threshold minimum weight between nodes in one component
         * */
        void connectedComponents(vector < set < unsigned  int  >> & components , unsigned  int threshold);
        /**
         * iteratively increases the weight threshold and each time finds the connected components which are assumed to be stable cores after doing multiple times of clustering
         * @param numOfClustering , number of clustering time
         * @param optimalClusterSize the optimal size of each cluster
         * @param numNodes number of ndoes in the graph (peired reads)
         * */
        vector < set < unsigned  int  >> getStableCore( size_t  numOfClustering , unsigned int optimalClusterSize ,  unsigned  int  numNodes);
        /**
         * prints nodes and edges in the graph
         **/
        void print( string name = "cytoNetwork1");
        
        /**
         * returns number of nodes in the graph. We initially added one to this value
         * 
         * */
        unsigned int getNumOfNodes(){
                return (numNodes-1);
        }
        /**
         * return the edge list of a source node
         * @param source the source node
         * @return the list of edges contains a target + weight
         * */
        vector<pair < unsigned  int  , char >> getEdgeList (unsigned  int source) {
                return adj[source];
        }
        /**
         * append the new list of edges to the current list
         * @param source the source node
         * @edgeList the list of edges should be appended to the current list
         * 
         * */
        void appendEdgeList (unsigned  int source, vector<pair < unsigned  int  , char >> edgeList) {
                adj[source].reserve(adj[source].size() + distance(edgeList.begin(),edgeList.end()));
                adj[source].insert (adj[source].end(), edgeList.begin(),edgeList.end());
        }
        
        
        
}; 
/**
 * This class builds the graph based on multiple clustering files. At the end, edges between nodes indicates they have been clustered at least once in the * same class. The weight of that edge is frequency of times that those two objects have been clustered together.
 * 
 * After building the Gragh , extractClusters finds the disconnected components in the graph whcih are the objects that have been clustered always together 
 * the procedure iteratively increase the weight threshold to finde the disconnected components (edges lower than threashold are ignored). Then the core of    * network can be find at higher value of threashold. For examle at weight = 8, if we find a disconnected component all the members of this component have * been clustered together at least 8 different times
 * */
class graphHandler {


private:
        undGraph G;                                     // The final graph which is built based on the thread graphs
        undGraph simNetwork;
        size_t numThreads;                              // Number of threads
        unsigned int numOfNodes;                        // Number of nodes in the graph
        /**
         * reads the clustering file an fill the dictionary data structure
         * @param fileName the clustering file name
         * @param clusterReadMap it keeps for each read to wich cluster belongs
         * 
         **/
        size_t getNumOfClusters( string  fileName);
        
        void  readAdjFromFile( string  fileName,  ClusterReadAss& clusterReadMap );
        /**
         *  @param threadID the ID of the current worker thread
         *  @param clusterReadMap A copy, it keeps for each read to wich cluster belongs
         *  @param threadGraph a graph wich belongs to this thread and stores part of edges accordign to his threadID
         * */
        void  addEdgeThread(size_t threadID ,const ClusterReadAss &clusterReadMap, undGraph &threadGraph );
        /**
         * builds the graph in parallel based on different clustering 
         * @param  fileNames clustering files 
         * */
        void  buildGraph(  vector <string> fileNames);
        
public :
        /**
         * finds the clusters in the range around of optimalClusterSize , we increase the range of accepted cluster size gradually
         * @param outputClusteringFile a file keeps the results , like final clustering 
         * @param iterNum number of iteration or number of different clustering
         * @param optimalClusterSize the optimal size of each cluster
         **/
        void  extractClusters ( string outputClusteringFile, size_t iterNum , unsigned int optimalClusterSize  );
        
        int getOptimalClusterSize(double coverage, string staFileName);
        
        void rescueUnclusterReads( std::map <unsigned int , unsigned int >   &readCluster);
        /**
         * constructor 
         * @param numOfNodes Number of nodes in the graph
         * @param fileNames list of different clustering files 
         * */
        graphHandler (unsigned int numOfNodes , vector <string> fileNames );
};
