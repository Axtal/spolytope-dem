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

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    DEM::Domain dom;
    dom.Alpha  = 3.0;
    dom.CamPos = Vec3_t(0.0,15.0,0.0); 

    dom.AddCylinder(-1, /*X0*/Vec3_t(-10.0,0.0,-2.0), /*R0*/0.5, /*X1*/Vec3_t(0.0,0.0,0.1), /*R1*/3.0, /*SR*/0.2, /*rho*/1.0);
    dom.GetParticle(-1)->FixVeloc();

    dom.AddTetra(-2,/*X*/Vec3_t(5.0,0.0,0.0),/*SR*/0.2,/*L*/2.0,/*rho*/3.0,0.0,&OrthoSys::e1);
    dom.GetParticle(-2)->v = -1.0,0.0,0.0;
    dom.GetParticle(-2)->w =  0.0,0.0,0.0;

    // Define the interaction constants
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu",1.0e5,1.0e5,8.0,8.0,0.0);
    B.Set(-2,"Kn Kt Gn Gt Mu",1.0e5,1.0e5,8.0,8.0,0.0);
    dom.SetProps(B);

    dom.Solve(/*Tf*/6.0e-3,/*dt*/3.0e-3,/*dtOut*/0.3,NULL,NULL,/*filekey*/"cone",/*RenderVideo?*/true);
    return 0;
}
MECHSYS_CATCH
