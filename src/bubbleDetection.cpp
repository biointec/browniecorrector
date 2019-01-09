/************************************************************************************
*    Copyright (C) 2014-2018 Jan Fostier (jan.fostier@ugent.be)                     *    
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

#include "graph.h"
#include "graphaln.h"
#include "settings.h"
#include "correctgraph.h"
#include "alignment.h"
#include <queue>
#include <map>


using namespace std;

bool DBGraph::handleParallelPaths(const vector<NodeID>& pathA,
                                  const vector<NodeID>& pathB,
                                  double covCutoff, size_t maxMargLength)
{
        //Check if two path are really parallel paths that needs to be similar and equal in length
        /*NWAligner ali(2, 1, -1, -3);
        string pathAstr = getPathSeq(pathA);
        string pathBstr = getPathSeq(pathB);
        if (std::abs( int( pathAstr.length()- pathBstr.length()) )> 2 ){
                return false;
        }
        if ( pathAstr.length() >settings.getK() && ali.align(pathAstr,pathBstr).score<((int)min( pathAstr.length(),pathBstr.length() ) / 3)){
                return false;
        }*/


        if (removeErroneousNodes(pathA, pathB, covCutoff, maxMargLength))
                return true;

        return( removeErroneousArcs (pathA, pathB, covCutoff, maxMargLength));


}
bool DBGraph::removeErroneousNodes(const vector<NodeID>& pathA,
                                   const vector<NodeID>& pathB,
                                   double covCutoff, size_t maxMargLength){
        size_t firstA, lastA;
        getUniquePath(pathA, firstA, lastA);
        vector<NodeID> subPathA;
        if (firstA <= lastA){
                subPathA = vector<NodeID>(pathA.begin() + firstA,
                                          pathA.begin() + lastA + 1);
        }
        size_t firstB, lastB;
        getUniquePath(pathB, firstB, lastB);
        vector<NodeID> subPathB;
        if (firstB <= lastB)
                subPathB = vector<NodeID>(pathB.begin() + firstB, pathB.begin() + lastB + 1);

        if (subPathA.empty() || subPathB.empty() )
                return false;

        double covSubPathA = getPathAvgKmerCov(subPathA);
        double covSubPathB =  getPathAvgKmerCov(subPathB);

        const vector<NodeID>& lowCovPath = (covSubPathA < covSubPathB) ?
        subPathA : subPathB;
        const double& lowCov = (covSubPathA < covSubPathB) ?
        covSubPathA : covSubPathB;
        

        bool remove = true;
        /*if (lowCov == covSubPathA){
                if (covSubPathB < covSubPathA * 2)
                        remove = false;
        }else{
                if (covSubPathA < covSubPathB * 2)
                        remove = false;
        }*/

        if (remove && lowCov <= covCutoff ) {
                //removePath(lowCovPath);
                flagPath(lowCovPath);
                return true;
        }
        return false;
}
bool DBGraph::removeErroneousArcs(const vector<NodeID>& pathA,
                                  const vector<NodeID>& pathB,
                                  double covCutoff, size_t maxMargLength){

        size_t firstA, lastA;
        getUniquePath(pathA, firstA, lastA);
        size_t firstB, lastB;
        getUniquePath(pathB, firstB, lastB);


        // Remove final arc?
        double covFinalArcA = ((pathA.size() >= 2) && (lastA == pathA.size() - 2)) ?
        getSSNode(pathA[lastA]).getRightArc(pathA[lastA+1])->getCoverage() : covCutoff + 1;
        double covFinalArcB = ((pathB.size() >= 2) && (lastB == pathB.size() - 2)) ?
        getSSNode(pathB[lastB]).getRightArc(pathB[lastB+1])->getCoverage() : covCutoff + 1;
        if (min(covFinalArcA, covFinalArcB) <= covCutoff) {
                if (covFinalArcA < covFinalArcB){
                        removeArc(pathA[lastA], pathA[lastA+1]);
                        return true;
                        /*if (covFinalArcA *2 < covFinalArcB){
                                flagArc(pathA[lastA], pathA[lastA+1]);
                                return true;
                        }*/
                }
                else{
                        removeArc(pathB[lastB], pathB[lastB+1]);
                        return true;
                        /*if (covFinalArcB *2 <covFinalArcA ){
                                flagArc(pathB[lastB], pathB[lastB+1]);
                                return true;
                        }*/
                }

        }

        // Remove first arc?
        double covFirstArcA = ((pathA.size() >= 2) && (firstA == 1)) ?
        getSSNode(pathA[0]).getRightArc(pathA[1])->getCoverage() : covCutoff + 1;
        double covFirstArcB = ((pathB.size() >= 2) && (firstB == 1)) ?
        getSSNode(pathB[0]).getRightArc(pathB[1])->getCoverage() : covCutoff + 1;
        if (min(covFirstArcA, covFirstArcB) <= covCutoff) {
                if (covFirstArcA < covFirstArcB){
                        removeArc(pathA[0], pathA[1]);
                        return true;
                        /*if (covFirstArcA*2 < covFirstArcB){
                                flagArc(pathA[0], pathA[1]);
                                return true;
                        }*/
                }
                else{
                        removeArc(pathB[0], pathB[1]);
                        return true;
                        /*if (covFirstArcB*2 <covFirstArcA ){
                                flagArc(pathB[0], pathB[1]);
                                return true;
                        }*/

                }

        }

        return false;
}


bool DBGraph::bubbleDetection(NodeID srcID, vector<NodeID>& visited,
                              vector<NodeID>& prevNode,
                              vector<NodeID>& nodeColor,
                              double covCutoff,
                              size_t maxMargLength, size_t maxNodesVisited)
{
        if (getSSNode(srcID).getNumRightArcs() < 2)
                return false;
        priority_queue<PathDFS, vector<PathDFS>, PathDFSComp> heap;
        heap.push(PathDFS(srcID, 0, 0));

        bool returnValue = false;

        while(!heap.empty()) {
                PathDFS currTop = heap.top();
                heap.pop();
                NodeID currID = currTop.nodeID;
                size_t currLength = currTop.length;
                SSNode curr = getSSNode(currID);
                if (!curr.isValid())
                        continue;
                for (ArcIt it = curr.rightBegin(); it != curr.rightEnd(); it++ ) {
                        NodeID nextID = it->getNodeID();
                        SSNode next = getSSNode(nextID);

                        // do we encounter a node previously encountered?
                        if ((prevNode[nextID + numNodes] != 0) || (nextID == srcID)) {
                                if (nodeColor[nextID + numNodes] == nodeColor[currID + numNodes])
                                        continue;

                                vector<NodeID> pathA = getPath(currID, prevNode);
                                pathA.push_back(nextID);
                                vector<NodeID> pathB = getPath(nextID, prevNode);

                                // if at least one node was deleted: get out of here!
                                if (handleParallelPaths(pathA, pathB, covCutoff, maxMargLength)) {
                                        returnValue = true;
                                        goto exitRoutine;
                                }
                        } else {
                                prevNode[nextID + numNodes] = currID;
                                nodeColor[nextID + numNodes] = (currID == srcID) ?
                                        nextID : nodeColor[currID + numNodes];
                                visited.push_back(nextID);

                                size_t nextLength = currLength + next.getMarginalLength();
                                if (nextLength > maxMargLength)
                                        continue;
                                if (visited.size() > maxNodesVisited)
                                        continue;
                                PathDFS nextTop(nextID, 0, nextLength);
                                heap.push(nextTop);
                        }
                }
        }

        // label definition to break out of nested loops
        exitRoutine:

        for (auto it : visited) {
                prevNode[it + numNodes] = 0;
                nodeColor[it + numNodes] = 0;
        }
        visited.clear();

        return returnValue;
}

void DBGraph::bubbleDetectionThread(size_t threadID, ParGraph& wlb,
                                    double covCutoff, size_t maxMargLength)
{
        vector<NodeID> visited;
        vector<NodeID> prevNode(2*numNodes+1, 0);
        vector<NodeID> nodeColor(2*numNodes+1, 0);

        while (true) {
                size_t firstNode, numNodes;
                wlb.getNodeChunk(firstNode, numNodes);

                if (numNodes == 0)
                        break;

                //cout << "Work from " << threadID << " " << firstNode << " to " << firstNode + numNodes << endl;
                for (size_t id = firstNode; id < firstNode + numNodes; id++) {
                        SSNode node = getSSNode(id);
                        if (!node.isValid())
                                continue;
                        // handle the positive node
                        bubbleDetection(id, visited, prevNode, nodeColor,
                                        covCutoff, maxMargLength,
                                        settings.getBubbleDFSNodeLimit());

                        // handle the negative node
                        bubbleDetection(-id, visited, prevNode, nodeColor,
                                        covCutoff, maxMargLength,
                                        settings.getBubbleDFSNodeLimit());
                }
        }
}

bool DBGraph::bubbleDetection(double covCutoff, size_t maxMargLength)
{
        const unsigned int& numThreads = settings.getNumThreads();
        ParGraph wlb(numNodes, settings.getThreadBubbleWorkSize());

        for (NodeID id = 1; id <= numNodes; id++)
                getSSNode(id).setFlag(false);

        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&DBGraph::bubbleDetectionThread, this,
                                          i, ref(wlb), covCutoff, maxMargLength);

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        cout << "\tProcessing graph (100%) " << endl;
        bool returnValue = false; size_t numNodesRemoved = 0;

        for (NodeID id = 1; id <= numNodes; id++) {
                if (getSSNode(id).getFlag()) {
                        removeNode(id);
                        //cout << "node " <<id << " detected as a bubble and removed" <<endl;

                        returnValue = true;
                        numNodesRemoved++;
                }
        }
        cout << "\tRemoved " << numNodesRemoved << " nodes" << endl;

        size_t numArcsRemoved = 0;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;

                vector<pair<NodeID, NodeID> > toRemove;

                for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                        if (node.getRightArc(it->getNodeID())->getFlag())
                                toRemove.push_back(pair<NodeID, NodeID>(id, it->getNodeID()));
                }

                for (auto it : toRemove) {
                        removeArc(it.first, it.second);
                        returnValue = true;
                        numArcsRemoved++;
                }
        }
        cout << "\tRemoved " << numArcsRemoved << " arcs" << endl;

        return numNodesRemoved >0;
}
