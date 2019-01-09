
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
#include "undGraph.h"
#include "calSim.h"
#include "extractReads.h"
#include "kmerObserver.h"
using namespace std;





int main(int argc, char ** args)
{
        if (argc <2){
                cout << "Wrong usage" <<endl;
                exit(0);
                
        }

        enum command {help, extractReadsCom, extractKmersCom, calSimNetworkCom, detectStableCommunityCom ,extractClustersCom , kmerObserverCom};
      
        command c;
        string arg(args[1]);
        if ((arg == "-h") || (arg == "--help")) {
                c = help;
        }else if (arg == "detectStableCommunity"){
                c = detectStableCommunityCom;
        }else if (arg == "calSim"){
                c = calSimNetworkCom;
        }else if (arg == "extractReads"){
                c = extractReadsCom;
        }else if (arg == "kmerObserver"){
             c = kmerObserverCom;       
        }else if (arg =="extractKmers"){
                c = extractKmersCom;
        }
        else{
                c = help;
        }
        switch (c)
        {
                case  detectStableCommunityCom:{
                        if (argc <8){
                                cout << "Correct usage:"<<endl;
                                cout << "./connectedComponents detectStableCommunity outputClusteringFile coverage staFile numOfNodes inputClusteringFile1 inputClusteringFile2 ... " << endl;
                        }
                        string outFileName = "";
                        string readFileName = "";
                        string stafileName = "";
                        double  coverage = 0 ;
                        
                        int numOfNodes = 0;
                        vector <string> fileNames;
                        for (size_t i = 2; i < argc; i++){
                                if (i == 2){
                                        outFileName = args[i];
                                        cout << "out FileName :\t" << outFileName <<endl; 

                                } else if (i == 3){
                                        coverage = stof (args[i]) ; 
                                        cout << "The coverage :\t" << coverage <<endl;
                                } else if (i == 4){
                                        stafileName =  args[i] ; 
                                        cout << "The sta file name :\t" << stafileName <<endl;
                                } else if (i == 5){
                                        numOfNodes = stoi (args[i]) ;
                                        cout << "Number Of Nodes :\t" << numOfNodes <<endl;
                                }
                                else if (i > 5){
                                        cout << endl;
                                        while (i < argc){
                                                fileNames.push_back(args[i]);
                                                cout << args[i]  << "   ";
                                                i = i +1;
                                        }
                                        cout <<endl;
                                        
                                }
                        }
                         Util::startChrono();
                         //findStableCommunityInRangeFromFile( numOfNodes , fileNames, optimalClusterSize , outFileName );
                          
                          graphHandler gh(numOfNodes ,fileNames );
                          int optimalClusterSize = gh.getOptimalClusterSize(coverage, stafileName);
                          gh.extractClusters(outFileName,fileNames.size(),optimalClusterSize);
                          cout << "Extracting clusters finished ("
                                << Util::stopChronoStr() << ")" << endl;
                        break;
                }
                case extractClustersCom :{
                        cout << "will be implemented soon" << endl;
                        break;
                };
                case help :{
                        cout << "./connectedComponents detectStableCommunity outputClusteringFile optimalClusterSize numOfNodes inputClusteringFile1 inputClusteringFile2 ...\n ";

                        break ;
                        
                };
                case calSimNetworkCom :
                {
                        if (argc < 6) {
                                cout << "ERROR: Wrong number of arguments!\n";
                                cout << "calSim brownie.15/2.corr.ncf network.dat brownie.15/nodes.stage3 15 readNodeAss.txt 30\n";
                                return 0;
                        }
                        string readFileName= args[2];
                        string ncfFile = args[3];
                        string networkFile = args[4];
                        string nodeFile = args[5];
                        size_t kmerSize = atoi( args[6]);
                        string readNodeAssFile = args[7];
                        double initialCov = atoi (args[8]);
                        string staFile =  (args[9]);
                        calSim cal (readFileName,ncfFile, nodeFile, kmerSize, readNodeAssFile , initialCov, staFile);
                        
                        cal.calSimInParallel(networkFile);
                        //findMutualNodesBetweenReads (readNodeSetRefined, networkFile, nodeLenVec , kmerSize);
                        break;
                        
                }
                case extractReadsCom :{
                    if (argc < 5) {
                                cout << "ERROR: Wrong number of arguments!\n";
                                cout << "extractReads kmerFile(fasta format) readFile \n";
                                return 0;
                        }
                        string kmerFile = args[2];
                        string readFile = args[3];
                        size_t kmerSize =atoi( args[4]);
                        string tempDir = args[5];
                        Util::startChrono();
                        extractReadsHandler er (kmerFile, readFile, kmerSize, tempDir);
                       
                           cout << "Separation of low comlex and normal reads finished in ("
                                << Util::stopChronoStr() << ")" << endl;
                        break;
                        
                }
                
        
                case kmerObserverCom : {
                        
                        if (argc < 5) {
                                cout << "ERROR: Wrong number of arguments!\n";
                                cout << "kmerObserver readFile kmerFile(fasta format) kmersize tempDir \n";
                                return 0;
                        }
                        string readFile = args[2];
                        string kmerFile = args[3];
                        size_t kmerSize = atoi( args[4]);
                        string tempDir = args[5];
                        extractKmersHandler er (kmerFile, readFile , kmerSize, tempDir);
                        er.extractKmerFreqFromFile(readFile);

                        
                }
               /* case clusteringCom :{
                        if (argc < 6) {
                                cout << "ERROR: Wrong number of arguments!\n";
                                cout << "clustering networkFile(input) clusteringFile(output) numOfNodes maxCutOffThreashold optimalClusterSize \n";
                                return 0;
                        }
                        string networkFile = args[2];
                        string outFileName = args[3];
                        unsigned int numOfNodes =  atoi( args[4]);                         
                        float maxCutOffThreashold = stof (args[5]) ;
                        unsigned int optimalClusterSize = atoi (args[6]) ;
                        undGraph g (numOfNodes , networkFile);
                        g.hierarchicalClustering(numOfNodes,maxCutOffThreashold,optimalClusterSize,outFileName );
                }  */
                
        }
        
        cout << "Finished successfully !!" <<endl;
        return 0;
        
        
}
