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

#ifndef KMERNPP_H
#define KMERNPP_H

#include "global.h"

// ============================================================================
// NODE POSITION PAIR CLASS
// ============================================================================

class NodePosPair : public std::pair<NodeID, NodePosition>
{
public:
        /**
         * Default constructor
         */
        NodePosPair() :
                std::pair<NodeID, NodePosition>(0, 0) {}

        /**
         * Default constructor
         * @param id Node identifier
         * @param pos Position identifier
         */
        NodePosPair(NodeID id, NodePosition pos) :
                std::pair<NodeID, NodePosition>(id, pos) {}

        /**
         * Get the node identifier
         * @return Node identifier
         */
        NodeID getNodeID() const {
                return first;
        }

        /**
         * Check whether the object points to valid position
         * @return True of false
         */
        bool isValid() const {
                return first != 0;
        }

        /**
         * Get the offset position
         * @return Offset position
         */
        NodePosition getPosition() const {
                return second;
        }
};

#endif
