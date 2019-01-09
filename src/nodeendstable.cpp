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

#include "nodeendstable.h"

bool NodeEndTable::insert(const Kmer& kmer, NodeID nodeID)
{
        assert(nodeID != 0);    // forbidden value

        // chose a representative kmer
        Kmer kmerRC = kmer.getReverseComplement();
        bool reverse = (doubleStranded) && (kmerRC < kmer);
        const Kmer &repKmer = (reverse) ? kmerRC : kmer;

        if (reverse)
                nodeID = -nodeID;

        return table.insert(Value(repKmer, NodeEndMD(nodeID))).second;
}

NodeEndRef NodeEndTable::find(const Kmer &kmer) const
{
        // chose a representative kmer
        Kmer kmerRC = kmer.getReverseComplement();
        bool reverse = (doubleStranded) && (kmerRC < kmer);
        const Kmer &repKmer = (reverse) ? kmerRC : kmer;

        return NodeEndRef(table.find(repKmer), reverse);
}
