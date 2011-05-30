/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    DEM::Domain dom;
    dom.CamPos = Vec3_t(0.0,10.0,0.0); 

    dom.AddCylinder (-4, /*X0*/Vec3_t(-10.0,0.0,0.0), /*R0*/0.01, /*X1*/Vec3_t(-3.0,0.0,0.1), /*R1*/3.0, 0.2, 1.0);
    dom.GetParticle(-4)->FixVeloc();

    dom.AddSphere(-1,Vec3_t(5.0,0.0,0.0),1.0,1.0);
    dom.GetParticle(-1)->v = -3.0,0.0,0.0;

    // properties of particles prior the triaxial test
    double dt = 1.0e-3;

    dom.Solve(10.0,dt,0.1,NULL,NULL,"cone",true);
    return 0;
}
MECHSYS_CATCH
