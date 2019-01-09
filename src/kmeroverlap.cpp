/************************************************************************************
*    Copyright (C) 2014-2016 Jan Fostier (jan.fostier@ugent.be)                     *    
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

#include "kmeroverlap.h"

using namespace std;

// ============================================================================
// KMER OVERLAP CLASS
// ============================================================================

const unsigned char KmerOverlap::leftMask[4] = {128, 64, 32, 16};
const unsigned char KmerOverlap::rightMask[4] = {8, 4, 2, 1};
const unsigned char KmerOverlap::cLeftMask = 240;
const unsigned char KmerOverlap::cRightMask = 15;
