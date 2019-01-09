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
#include <string>

TEST(Scaffold, CircAPrioriSnap)
{
       /* Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(2, 3), GaussVal(-500, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());*/
}

TEST(Scaffold, CircTest)
{
      /*  Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(3), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(200, 10), scaffold.getNodeInfo(3).first);*/
}

TEST(Scaffold, ReverseComplPresent)
{
   /*     Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.addObservation(NodePair(1, 4), GaussVal(300, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(2, -3), GaussVal(100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(2), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(300, 10), scaffold.getNodeInfo(4).first);*/
}

TEST(Scaffold, SnapAPrioriLinkAtRootNode)
{
     /*   Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(200, 10));
        scaffold.calculateScaffold();

        scaffold.addObservation(NodePair(1, 3), GaussVal(1, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);*/
}

TEST(Scaffold, SnapPosterioriLinkAtRootNode)
{
     /*   Scaffold scaffold(1);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(400, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);*/
}

TEST(Scaffold, Linear)
{
    /*    Scaffold scaffold(1);

        EXPECT_EQ(size_t(1), scaffold.getScaffold().size());
        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(1).first);

        scaffold.addObservation(NodePair(1, 2), GaussVal(100, 10));
        scaffold.addObservation(NodePair(1, 3), GaussVal(-100, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(scaffold.getScaffold().size(), size_t(3));
        EXPECT_EQ(GaussVal(100, 10), scaffold.getNodeInfo(2).first);
        EXPECT_EQ(GaussVal(-100, 10), scaffold.getNodeInfo(3).first);
        EXPECT_EQ(scaffold.getNodeInfo(1).first, GaussVal(0, 0));

        scaffold.addObservation(NodePair(3, 2), GaussVal(200, 10));
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(100, 1.0/sqrt(0.015)), scaffold.getNodeInfo(2).first);
        EXPECT_EQ(GaussVal(-100, 1.0/sqrt(0.015)), scaffold.getNodeInfo(3).first);

        scaffold.setRootID(3);
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(3).first);

        scaffold.setRootID(2);
        scaffold.calculateScaffold();

        EXPECT_EQ(GaussVal(0, 0), scaffold.getNodeInfo(2).first);*/
}


