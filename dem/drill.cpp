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
    dom.CamPos = Vec3_t(0.0,20.0,3.0); 
    dom.GenSpheres  (-1, 10, 10, 3.0, "HCP", 1000, 0.9,1.0);
    dom.GenOpenBox (-2,15.0,15.0,10.0,0.1,1.0);
    dom.GetParticle(-2)->FixVeloc();
    dom.GetParticle(-3)->FixVeloc();
    dom.GetParticle(-4)->FixVeloc();
    dom.GetParticle(-5)->FixVeloc();
    dom.GetParticle(-6)->FixVeloc();

    Vec3_t pos(0.0,0.0,10.0);
    dom.AddDrill (-7,pos,0.1,5.0,10.0,3.0);
    dom.GetParticle(-7)->FixVeloc();

    dom.WriteBPY("drill");

    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Initialize(i);
        dom.Particles[i]->Ff = 0.0,0.0,-dom.Particles[i]->Props.m*9.8;
    }

    dom.Solve(10.0,1.0e-3,0.1,NULL,NULL,"drilla",true);

    dom.GetParticle(-7)->vzf=false;
    
    dom.Solve(20.0,1.0e-3,0.1,NULL,NULL,"drillb",true);
    return 0;
}
MECHSYS_CATCH
