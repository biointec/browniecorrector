/************************************************************************************
*    Copyright (C) 2015 Jan Fostier (jan.fostier@ugent.be)                          *    
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
#include <gtest/gtest.h>
#include <cstdio>
#include "tkmer.h"

using namespace std;

TEST(nucleotide, basicTest)
{
        EXPECT_EQ(Nucleotide::getComplement('A'), 'T');
        EXPECT_EQ(Nucleotide::getComplement('T'), 'A');
        EXPECT_EQ(Nucleotide::getComplement('G'), 'C');
        EXPECT_EQ(Nucleotide::getComplement('C'), 'G');
}

TEST(nucleotide, stringTest)
{
        string str("AGCTTGCT");

        Nucleotide::reverse(str);
        EXPECT_EQ(str == "TCGTTCGA", true);
        Nucleotide::complement(str);
        EXPECT_EQ(str == "AGCAAGCT", true);
        Nucleotide::revCompl(str);
        EXPECT_EQ(str == "AGCTTGCT", true);
}

TEST(nucleotide, NWTest)
{
        const int mSize = 100;

        string seq1 = ("AGCTTGTAGCTTGCTAGCTTGCTAGCTTGCT");
        string seq2 = ("AGCAAGCTAGCTTGCTAGCCTAGCTTGCT");
        //              AGCTTG-TAGCTTGCTAGCTTGCTAGCTTGCT
        //              AGCAAGCTAGCTTGCTAGC---CTAGCTTGCT

        // allocate alignment matrix
        int **A = new int*[mSize];
        for (int i = 0; i < mSize; i++)
                A[i] = new int[mSize];

        for (int i = 0; i < mSize; i++) {
                A[i][0] = i*indelScore;
                A[0][i] = i*indelScore;
        }

        Nucleotide::getNWAlignMatrix(seq1, seq2, A);

        EXPECT_EQ(A[seq1.size()][seq2.size()], 26);

        int maxScore = Nucleotide::maxNWScore(seq1);
        EXPECT_EQ((size_t)maxScore, seq1.size());

        string align1, align2;
        Nucleotide::getNWAlignment(seq1, seq2, A, align1, align2);

        EXPECT_EQ(align1, string("AGCTTG-TAGCTTGCTAGCTTGCTAGCTTGCT"));
        EXPECT_EQ(align2, string("AGCAAGCTAGCTTGCTAGC---CTAGCTTGCT"));

        for (int i = 0; i < mSize; i++)
                delete [] A[i];
        delete [] A;

}
