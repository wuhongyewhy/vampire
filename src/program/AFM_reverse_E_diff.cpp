//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify 
//  it under the terms of the GNU General Public License as published by 
//  the Free Software Foundation; either version 2 of the License, or 
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but 
//  WITHOUT ANY WARRANTY; without even the implied warranty of 
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License 
//  along with this program; if not, write to the Free Software Foundation, 
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
///
/// @file
/// @brief calculate low energy state for material named "AFMsc" 
///
/// @details make magnetization for material named "AFMsc" reverse,
///  and calculate energy different. Chose a low energy state.
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section make magnetization for material named "AFMsc" reverse,
///  and calculate energy different. Chose a low energy state.
/// @author  Wu Hong-ye, wuhongyewhy@hotmail.com
/// @version 1.0
/// @date    19/06/2014
/// @internal
///	Created:		19/06/2014
///	Revision:	 
///=====================================================================================
///

// Standard Libraries
#include <iostream>
#include <sstream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace program{
	
void AFM_low_E_state(){
	// check calling of routine if error checking is activated
	if(err::check==true){std::cout << "program::AFM_low_E_state has been called" << std::endl;}

    bool last_state = stats::calculate_energy;
    stats::calculate_energy = true;
    stats::mag_m();

    std::stringstream stream;
    stats::output_energy(stream, stats::all, stats::total);

    double energy_old;
    stream >> energy_old;

    cs::AFM_reverse();
    stats::mag_m();
    stats::calculate_energy = last_state;

    stream.clear();
    stats::output_energy(stream, stats::all, stats::total);
    double energy_new;
    stream >> energy_new;
    stream.clear();

    terminaltextcolor(GREEN);
    std::cout << "Old AFM state energy: " << energy_old << "\n";
    std::cout << "New AFM state energy: " << energy_new << "\n";

    if (energy_old < energy_new) {
        cs::AFM_reverse();
        std::cout << "Chosing ";
        terminaltextcolor(RED);
        std::cout << "OLD";
        terminaltextcolor(GREEN);
        std::cout << " state\n";
    } else {
        std::cout << "Chosing ";
        terminaltextcolor(RED);
        std::cout << "NEW";
        terminaltextcolor(GREEN);
        std::cout << " state\n";
    }
    terminaltextcolor(WHITE);
	
	return;
}

}//end of namespace program

