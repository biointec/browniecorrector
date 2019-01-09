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
#include "util.h"

TEST(poissonPDF, poissonPDFTest)
{
        EXPECT_DOUBLE_EQ(0.12511003572113349, Util::poissonPDF(10, 10));
        EXPECT_DOUBLE_EQ(4.8646491820674864E-63, Util::poissonPDF(100, 10));
        EXPECT_DOUBLE_EQ(0, Util::poissonPDF(1000, 10));
        EXPECT_DOUBLE_EQ(0, Util::poissonPDF(10000, 10));
}
