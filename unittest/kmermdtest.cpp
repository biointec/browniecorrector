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
#include "kmeroverlap.h"
#include "tkmer.h"

TEST(kmerMD, kmerMDTest)
{
        KmerOverlap md;
        char nucleotide;

        EXPECT_EQ(md.hasLeftOverlap('A'), false);
        md.markLeftOverlap('A');
        EXPECT_EQ(md.hasUniqueLeftOverlap(nucleotide), true);
        EXPECT_EQ(nucleotide, 'A');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), false);
        EXPECT_EQ(md.hasLeftOverlap('A'), true);
        md.markRightOverlap('C');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), true);
        EXPECT_EQ(nucleotide, 'C');
        md.markRightOverlap('G');
        EXPECT_EQ(md.hasUniqueRightOverlap(nucleotide), false);
        EXPECT_EQ(md.hasRightOverlap('C'), true);
        EXPECT_EQ(md.hasLeftOverlap('G'), false);
}
