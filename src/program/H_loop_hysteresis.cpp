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
/// @brief Contains the H_loop_hysteresis program
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section info File Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///

// Standard Libraries
#include <iostream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

namespace program{

/// @brief Function to calculate magnetisation over a time series
///
/// @section License
/// Use of this code, either in source or compiled form, is subject to license from the authors.
/// Copyright \htmlonly &copy \endhtmlonly Richard Evans, 2009-2010. All Rights Reserved.
///
/// @section Information
/// @author  Richard Evans, richard.evans@york.ac.uk
/// @version 1.0
/// @date    30/03/2011
///
/// @internal
///	Created:		30/03/2011
///	Revision:	--
///=====================================================================================
///
void H_loop_hysteresis(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::H_loop_hysteresis has been called" << std::endl;

	double temp=sim::temperature;
	
	// Set equilibration temperature and field
	sim::temperature=sim::Teq;
	
	// Equilibrate system
	while(sim::time<sim::equilibration_time){
		
		sim::integrate(sim::partial_time);
		
	}
	// Calculate magnetisation statistics
	stats::mag_m();
		
	// Output data
	vout::data();
    // atoms output
    if ((vout::output_atoms_config==true)){
        if (vout::output_atoms_file_counter==0) vout::atoms_coords();
        vout::atoms();
    }

	
    sim::Hmax = abs(sim::Hmax);
    sim::Hinc = abs(sim::Hinc);
	sim::temperature=temp;
    sim::H_applied = sim::Hmax;
	
    double H_delta;
	// Perform Time Series
    while (sim::H_applied > -sim::Hmax){

		// Integrate system
		sim::integrate(1);
		
        H_delta = sim::Hinc*(sim::time-sim::equilibration_time);
        sim::H_applied = sim::Hmax - H_delta;


        if (vout::output_atoms_config && 
            ((sim::time-sim::equilibration_time)%vout::output_atoms_config_rate==0)) {
            // Calculate magnetisation statistics
            stats::mag_m();

            // Output data
            vout::data();
            // atoms output
            if ((vout::output_atoms_config==true)){
                if (vout::output_atoms_file_counter==0) vout::atoms_coords();
                vout::atoms();
            }
        }

	}

    int timewall = sim::time;

    while (sim::H_applied < sim::Hmax){

        // Integrate system
        sim::integrate(1);

        H_delta = sim::Hinc*(sim::time-timewall);
        sim::H_applied = -sim::Hmax + H_delta;


        if (vout::output_atoms_config &&
            ((sim::time-sim::equilibration_time)%vout::output_atoms_config_rate==0)) {
            // Calculate magnetisation statistics
            stats::mag_m();

            // Output data
            vout::data();
            // atoms output
            if ((vout::output_atoms_config==true)){
                if (vout::output_atoms_file_counter==0) vout::atoms_coords();
                vout::atoms();
            }
        }

    }

}

}//end of namespace program
