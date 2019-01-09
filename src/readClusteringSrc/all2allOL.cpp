/************************************************************************************
*    Copyright (C) 2018 Jan Fostier (jan.fostier@ugent.be)                 *    
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


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>
#include <mutex>
#include <math.h>
#include <iomanip>

using namespace std;

// ============================================================================
// AUXILIARY ROUTINES
// ============================================================================

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

// ============================================================================
// WORK LOAD BALANCER
// ============================================================================

class WorkLoadBalancer
{
private:
        size_t matrixDim;               // matrix dimension
        size_t chunkSize;               // target chunk size
        size_t iCurr;                   // current i-index

        std::mutex mutex;               // mutex
        std::string message;            // message to print during progress

public:
        /**
         * Get a chunk of nodes (thread-safe)
         * @param iFirst First i-index to handle
         * @param iEnd End i-index to handle
         * @return Number of elements to handle
         */
        bool getChunk(size_t& iFirst, size_t& iEnd) {
                // lock the mutex (RAII)
                std::unique_lock<std::mutex> lock(mutex);

                // no more work: send termination signal
                if (iCurr >= matrixDim)
                        return false;

                // assign a workload
                iFirst = iCurr;
                iCurr = min<size_t>(iCurr + chunkSize, matrixDim);
                iEnd = iCurr;

                double perc = 100.0 * (double)iCurr / (double)matrixDim;
                cout << fixed << setprecision(1) << message << " (" << perc << "%)";
                if (iCurr >= matrixDim) {
                        cout << endl;
                } else {
                        cout << "\r";
                        cout.flush();
                }

                return true;
        }

        /**
         * Default constructor
         * @param matrixDim Matrix dimension
         * @param chunkSize Chunk size
         * @param message Message to print during progress
         */
        WorkLoadBalancer(size_t matrixDim, size_t chunkSize,
                         const std::string& message = "Processing") :
                matrixDim(matrixDim), chunkSize(chunkSize), iCurr(0),
                message(message) {}
};

// ============================================================================
// MATRIX CLASS
// ============================================================================

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

// ============================================================================
// READ PAIR CLASS
// ============================================================================

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
};

// ============================================================================
// MAIN ROUTINES
// ============================================================================

/**
 * Write program usage information to the standard output
 */
void printUsage()
{
        cout << "Usage: all2allOL input.fastqecd  pattern.txt\n\n";
        cout << "Report bugs to Jan Fostier <jan.fostier@ugent.be>" << endl;
}

/**
 * Read sequences from a FASTQ file
 * @param filename FASTQ file filename
 * @param sequences Vector of sequences (output)
 */
void readSequences(const string& filename, vector<ReadPair>& readPairs)
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

void readPattern(const string& filename, string& pattern)
{
        ifstream ifs(filename.c_str());
        if (!ifs)
                throw runtime_error("Could not open file: " + filename);

        while (ifs) {
                getline(ifs, pattern);
                if (pattern.empty())
                        continue;
                break;
        }

        if (pattern.empty())
                throw runtime_error("File: " + filename + " should contain a pattern");
}

/**
 * Perform global alignment of two sequences and print the alignment to stdout
 * @param X sequence one
 * @param Y sequence two
 */
int alignOverlap(const string& X, const string& Y)
{
        const int G = -3;
        const int M = 1;
        const int I = -1;

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

int computeScore(const ReadPair& p1, const ReadPair& p2)
{
        // forward forward alignment
        int FF1 = alignOverlap(p1.read1, p2.read1) + alignOverlap(p1.read2, p2.read2);
        int FF2 = alignOverlap(p1.read2, p2.read1);
        int FF3 = alignOverlap(p1.read1, p2.read2);

        int maxFF = max<int>(FF1, max<int>(FF2, FF3));

        ReadPair p1RC;
        p1RC.read1 = revComp(p1.read2);
        p1RC.read2 = revComp(p1.read1);

        // forward forward alignment
        int RF1 = alignOverlap(p1RC.read1, p2.read1) + alignOverlap(p1RC.read2, p2.read2);
        int RF2 = alignOverlap(p1RC.read2, p2.read1);
        int RF3 = alignOverlap(p1RC.read1, p2.read2);

        int maxRF = max<int>(RF1, max<int>(RF2, RF3));

        return max<int>(maxFF, maxRF);
}

void workerThread(WorkLoadBalancer& wlb, const vector<ReadPair>& readPairs,
                  ofstream& ofs, mutex& m)
{
        vector<int> scores;

        size_t iBegin, iEnd;
        while (wlb.getChunk(iBegin, iEnd)) {
                // compute the scores and buffer them
                for (size_t i = iBegin; i < iEnd; i++)
                        for (size_t j = i+1; j < readPairs.size(); j++)
                                scores.push_back(computeScore(readPairs[i], readPairs[j]));

                // lock the mutex (RAII) and write results to disk
                std::unique_lock<std::mutex> lock(m);

                vector<int>::iterator s = scores.begin();
                for (size_t i = iBegin; i < iEnd; i++)
                        for (size_t j = i+1; j < readPairs.size(); j++)
                                ofs << i+1 << "\t" << j+1 << "\t" << *s++ << "\n";

                scores.clear();
        }
}

int main(int argc, char** argv)
{
        // check for the correct number of arguments
        if (argc != 3) {
                printUsage();
                return EXIT_FAILURE;
        }

        // load the read pairs from disk into memory
        vector<ReadPair> readPairs;
        readSequences(argv[1], readPairs);
        cout << "Loaded " << readPairs.size() << " read pairs" << endl;

        // load the pattern into memory
        string pattern;
        readPattern(argv[2], pattern);
        cout << "Loaded pattern: " << pattern << endl;

        // identify the patterns in the read pairs
        for (auto pair: readPairs)
                pair.findPattern(pattern);

        cout << "Succesfully verified that all read pairs contain the pattern" << endl;

        // start the threads
        ofstream ofs("distances.txt");
        mutex outputMutex;

        unsigned int numThreads = thread::hardware_concurrency();
        cout << "Computing all pairwise alignments using " << numThreads << " threads" << endl;

        WorkLoadBalancer wlb(readPairs.size(), 1);
        vector<thread> workerThreads(numThreads);
        for (size_t i = 0; i < workerThreads.size(); i++)
                workerThreads[i] = thread(&workerThread, ref(wlb), cref(readPairs),
                                          ref(ofs), ref(outputMutex));

        // wait for worker threads to finish
        for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

        ofs.close();

        cout << "Exiting... bye!" << endl;

        return EXIT_SUCCESS;
}
