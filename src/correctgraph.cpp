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

// ============================================================================
// PRIVATE CORRECTGRAPH.CPP (STAGE 4 ROUTINES)
// ============================================================================

void DBGraph::removeNode(NodeID nodeID)
{
/*#ifdef DEBUG
        if (trueMult[abs(nodeID)] > 0)
                cout << "\tERROR removing node " << nodeID << endl;
#endif*/

        SSNode node = getSSNode(nodeID);
        for (ArcIt it = node.leftBegin(); it != node.leftEnd(); it++) {
                SSNode leftNode = getSSNode(it->getNodeID());
                if (leftNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = leftNode.deleteRightArc (node.getNodeID());
                assert(result);
        }

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++) {
                SSNode rightNode = getSSNode ( it->getNodeID() );
                if (rightNode.getNodeID()==-node.getNodeID())
                        continue;
                bool result = rightNode.deleteLeftArc(node.getNodeID());
                assert(result);
        }

        node.deleteAllRightArcs();
        node.deleteAllLeftArcs();
        node.invalidate();
        numValidNodes--;
}

void DBGraph::flagNode(NodeID nodeID)
{
        getSSNode(nodeID).setFlag(true);
}

void DBGraph::removeArc(NodeID leftID, NodeID rightID)
{
#ifdef DEBUG
        if (getSSNode(leftID).getRightArc(rightID)->getTrueArc()) {
                cout << "\tERROR detaching nodes " << leftID << " and " << rightID << endl;
                //exit(EXIT_SUCCESS);
        }
#endif

        getSSNode(leftID).deleteRightArc(rightID);
        getSSNode(rightID).deleteLeftArc(leftID);
}

void DBGraph::flagArc(NodeID leftID, NodeID rightID)
{
        getSSNode(leftID).getRightArc(rightID)->setFlag(true);
        getSSNode(rightID).getLeftArc(leftID)->setFlag(true);
}

vector<NodeID> DBGraph::getPath(NodeID dstID, const vector<NodeID>& prevNode) const
{
        vector<NodeID> path;
        path.push_back(dstID);

        while (true) {
                dstID = prevNode[dstID + numNodes];
                if (dstID == 0)
                        break;
                path.push_back(dstID);
        }

        reverse(path.begin(), path.end());
        return path;
}

void DBGraph::getUniquePath(const std::vector<NodeID>& path,
                            size_t& first, size_t& last) const
{
        // sanity check is necessary
        if (path.size() <= 2) {
                first = 1; last = 0;
                return;
        }

        // find the first node that can be deleted
        first = 1;
        for (size_t i = path.size() - 2; i >= 1; i--) {
                if (getSSNode(path[i]).getNumRightArcs() == 1)
                        continue;
                first = i+1;
                break;
        }

        // find the last node that can be deleted
        last = path.size() - 2;
        for (size_t i = 1; i <= path.size() - 2; i++) {
                if (getSSNode(path[i]).getNumLeftArcs() == 1)
                        continue;
                last = i-1;
                break;
        }
}
string DBGraph::getPathSeq(const vector<NodeID>& path){
        string seq = "";
        for (size_t i=1;i<path.size()-1;i++) {
                seq = seq + getSSNode(path[i]).substr(Kmer::getK()-1, getSSNode(path[i]).getMarginalLength());
        }
        return seq;

}
double DBGraph::getPathAvgKmerCov(const vector<NodeID>& path)
{
        double totKmer = 0, totLen = 0;
        for (auto it : path) {
                totKmer += getSSNode(it).getKmerCov();
                totLen += getSSNode(it).getMarginalLength();
        }

        return totKmer / totLen;
}

void DBGraph::removePath(const vector<NodeID>& pathA)
{
        for (auto it : pathA)
                removeNode(it);
}

void DBGraph::flagPath(const vector<NodeID>& pathA)
{
        for (NodeID nodeID : pathA)
                flagNode(nodeID);
}

bool DBGraph::flowCorrection(NodeID nodeID)
{
        SSNode node = getSSNode(nodeID);
        if (!node.isValid())
                return false;

        int expNodeMult = getExpMult(node.getAvgKmerCov());
        if (expNodeMult == 0)
                return false;

        //cout << "Multiplicity for node " << nodeID << ": " << expNodeMult << endl;

        int sumArcMult = 0;
        bool candidateRemoval = false;
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0) {
                        candidateRemoval = true;
                        expArcMult++;
                }
                sumArcMult += expArcMult;
        }

        //cout << "Sum of the right arc multiplicities: " << sumArcMult << endl;

        // we will not detach arcs in this step
        if (sumArcMult <= expNodeMult && !candidateRemoval)
                return false;

        //cout << "Sum of arcs is higher than expected multiplicity" << endl;

        // a) First assume that the topology is CORRECT
        double totCorrProb = getObsProbLog(node.getAvgKmerCov(), node.getMarginalLength(), sumArcMult);
        //cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << sumArcMult << ": " << totCorrProb << endl;

        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)    // bring this to one as we assume topology to be correct
                        expArcMult++;
                double arcProb = getObsProbLog(it->getCoverage(), 1, expArcMult);
                //cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totCorrProb += arcProb;
        }

        //cout << "TOTAL log prob assuming topology is correct: " << totCorrProb << endl;

        // b) Now assume that the topology is INCORRECT
        double totWrongProb = getObsProbLog(node.getAvgKmerCov(), node.getMarginalLength(), expNodeMult);
        //cout << "Log prob of node " << nodeID << " with coverage: " << node.getAvgKmerCov() << "  having multiplicity: " << expNodeMult << ": " << totWrongProb << endl;

        vector<NodeID> toDetach;
        for (ArcIt it = node.rightBegin(); it != node.rightEnd(); it++ ) {
                int expArcMult = getExpMult(it->getCoverage());
                if (expArcMult == 0)
                        toDetach.push_back(it->getNodeID());
                double arcProb = getObsProbLog(it->getCoverage(), 1, expArcMult);
                //cout << "Log prob of arc with coverage " << it->getCoverage() << " having multiplicity: " << expArcMult << ": " << arcProb << endl;
                totWrongProb += arcProb;
        }

        //cout << "TOTAL log prob assuming topology is WRONG: " << totWrongProb << endl;

        if (totWrongProb - totCorrProb < 5)
                return false;

        for (auto it : toDetach) {
                removeArc(nodeID, it);
                if ((getSSNode(it).getNumLeftArcs() == 0) &&
                    (getSSNode(it).getNumRightArcs() == 0)) {
                        getSSNode(it).invalidate();
                        numValidNodes--;
                }
        }
        return !toDetach.empty();
}

// ============================================================================
// PUBLIC CORRECTGRAPH.CPP (STAGE 4 ROUTINES)
// ============================================================================


void DBGraph::concatenateAroundNode(NodeID seedID, vector<NodeID>& nodeListv)
{
        nodeListv.clear();

        SSNode seed = getSSNode(seedID);
        if (!seed.isValid())
                return;

        deque<NodeID> nodeListq;
        nodeListq.push_back(seedID);
        seed.setFlag(true);

        // find linear paths to the right
        SSNode curr = seed;
        while (curr.getNumRightArcs() == 1) {
                NodeID rightID = curr.rightBegin()->getNodeID();
                SSNode right = getSSNode(rightID);
                // don't merge palindromic repeats / loops
                if (right.getFlag())
                        break;
                if (right.getNumLeftArcs() != 1)
                        break;
                nodeListq.push_back(rightID);
                right.setFlag(true);
                curr = right;
        }

        // find linear paths to the left
        curr = seed;
        while (curr.getNumLeftArcs() == 1) {
                NodeID leftID = curr.leftBegin()->getNodeID();
                SSNode left = getSSNode(leftID);
                // don't merge palindromic repeats / loops
                if (left.getFlag())
                        break;
                if (left.getNumRightArcs() != 1)
                        break;
                nodeListq.push_front(leftID);
                left.setFlag(true);
                curr = left;
        }

        // reset the flags to false
        for (const auto& it : nodeListq)
                getSSNode(it).setFlag(false);

        // if no linear path was found, continue
        if (nodeListq.size() == 1)
                return;

        // concatenate the path
        NodeID frontID = nodeListq.front();
        SSNode front = getSSNode(frontID);
        NodeID backID = nodeListq.back();
        SSNode back = getSSNode(backID);

        front.deleteAllRightArcs();
        front.inheritRightArcs(back);

        size_t newKmerCov = front.getKmerCov();
        size_t newReadStartCov = front.getReadStartCov();
        for (size_t i = 1; i < nodeListq.size(); i++) {
                newKmerCov += getSSNode(nodeListq[i]).getKmerCov();
                newReadStartCov += getSSNode(nodeListq[i]).getReadStartCov();
                getSSNode(nodeListq[i]).deleteAllLeftArcs();
                getSSNode(nodeListq[i]).deleteAllRightArcs();
                getSSNode(nodeListq[i]).invalidate();
        }

        front.setKmerCov(newKmerCov);
        front.setReadStartCov(newReadStartCov);

        copy(nodeListq.begin(), nodeListq.end(), std::back_inserter(nodeListv));

        string str;
        convertNodesToString(nodeListv, str);

        front.setSequence(str);
        numValidNodes -= nodeListq.size() - 1;
}

bool DBGraph::concatenateNodes()
{
        size_t numConcatenations = 0;
#ifdef DEBUG
        size_t numIncorrectConcatenations = 0;
#endif

        for (NodeID seedID = 1; seedID <= numNodes; seedID++) {
                vector<NodeID> concatenation;
                concatenateAroundNode(seedID, concatenation);

                if (!concatenation.empty())
                        numConcatenations += concatenation.size() - 1;

#ifdef DEBUG

                if (trueMult.empty())
                        continue;

                for (size_t i = 1; i < concatenation.size(); i++) {
                        NodeID lID = concatenation[i-1];
                        NodeID rID = concatenation[i];
                        SSNode left = getSSNode(lID);
                        SSNode right = getSSNode(rID);
                        size_t lMult = trueMult[abs(lID)];
                        size_t rMult = trueMult[abs(rID)];
                        if (lMult != rMult) {
                                cout << "\tConcatenating " << lID << " (cov = "
                                     << left.getAvgKmerCov() << ", " << lMult
                                     << ") and " << rID << " (cov = "
                                     << right.getAvgKmerCov() << ", " << rMult
                                     << ")" << endl;
                                numIncorrectConcatenations++;
                                trueMult[abs(lID)] = max(lMult, rMult);
                                trueMult[abs(rID)] = max(lMult, rMult);
                        }
                }
#endif
        }

        cout << "\tConcatenated " << numConcatenations << " nodes" << endl;

#ifdef DEBUG
        if (numIncorrectConcatenations > 0)
                cout << "\t" << "Number of incorrect connections: "
                     << numIncorrectConcatenations << endl;
#endif

        size_t countTotal = 0;
        for (NodeID seedID = 1; seedID <= numNodes; seedID++)
                if (getSSNode(seedID).isValid())
                        countTotal++;

        return (numConcatenations > 0);
}


bool DBGraph::flowCorrection()
{
        bool returnValue = false;
        for (NodeID id = -numNodes; id <= numNodes; id++) {
                if (id == 0)
                        continue;
                if (!getSSNode(id).isValid())
                        continue;
                if (getSSNode(id).getNumRightArcs() < 2)
                        continue;

                if (abs(id) % OUTPUT_FREQUENCY == 0) {
                        cout << "\tProcessing node " << id << "/" << numNodes << "\r";
                        cout.flush();
                }

                if (flowCorrection(id))
                        returnValue = true;
        }

        cout << "\tProcessing node " << numNodes << "/" << numNodes << endl;

        return returnValue;
}

