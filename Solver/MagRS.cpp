/* 
 * MagRS.cpp, part of the RotationShield program
 * Copyright (c) 2007-2014 Michael P. Mendenhall
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "MagRS.hh"

BField_Protocol* BField_Protocol::BFP = new BField_Protocol();
    
//-----------------------------------------------------

bool MagExtField::queryInteraction(void* ip) {

    if(ip != BField_Protocol::BFP) return false;

    if(BField_Protocol::BFP->M2) {
        BField_Protocol::BFP->M2B += fieldAtWithTransform2(BField_Protocol::BFP->x, *BField_Protocol::BFP->M2);
    } else if(BField_Protocol::BFP->M3) {
        BField_Protocol::BFP->B += fieldAtWithTransform3(BField_Protocol::BFP->x, *BField_Protocol::BFP->M3);
    } else {
        BField_Protocol::BFP->B += fieldAt(BField_Protocol::BFP->x);
    }
    return true;
}
