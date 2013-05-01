/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Torres                                     *
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
    dom.CamPos = Vec3_t(0.0,30.0,6.0); 
    dom.GenSpheres  (-1, 10, 10, 3.0, "HCP", 1000, 0.9,1.0);
    dom.Center();
    dom.AddCylinder (-2, /*X0*/Vec3_t(0.0,0.0,-6.0), /*R0*/10.0, /*X1*/Vec3_t(0.0,0.0,4.0), /*R1*/10.0, 0.1, 3.0);
    dom.AddPlane    (-3, Vec3_t(0.0,0.0,-6.0), 0.1, 20.0, 20.0, 1.0);
    dom.GetParticle(-2)->FixVeloc();
    dom.GetParticle(-3)->FixVeloc();

    //dom.AddCylinder (-4, /*X0*/Vec3_t(0.0,0.0,11.0), /*R0*/3.0, /*X1*/Vec3_t(0.0,0.0,6.0), /*R1*/0.001, 0.2, 1.0);
    //dom.GetParticle(-4)->FixVeloc();
    //dom.GetParticle(-4)->vzf=false;

    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Initialize(i);
        dom.Particles[i]->Ff = 0.0,0.0,-dom.Particles[i]->Props.m*9.8;
    }

    // properties of particles prior the triaxial test
    double dt = 1.0e-3;
    double Kn = 1.0e5;
    double Kt = 5.0e4;
    double Gn = 16.0;
    double Gt = 8.0;
    double Mu = 0.4;
    Dict B;
    B.Set(-1,"Kn Kt Mu Gn Gt",Kn,Kt,Mu,Gn,Gt);
    B.Set(-2,"Kn Kt Mu Gn Gt",Kn,Kt,Mu,Gn,Gt);
    B.Set(-3,"Kn Kt Mu Gn Gt",Kn,Kt,Mu,Gn,Gt);
    B.Set(-4,"Kn Kt Mu Gn Gt",Kn,Kt,Mu,Gn,Gt);
    dom.SetProps(B);

    dom.Solve(10.0,dt,0.1,NULL,NULL,"drilla",true);
    
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    dom.AddCylinder(-4, /*X0*/Vec3_t(0.0,0.0,Xmax(2)+3.0+0.2), /*R0*/2.0, /*X1*/Vec3_t(0.0,0.0,Xmax(2)+0.2), /*R1*/0.001, 0.2, 3.0);
    dom.GetParticle(-4)->Initialize(dom.Particles.Size()-1);
    dom.GetParticle(-4)->InitializeVelocity(dt);
    dom.GetParticle(-4)->Ff = 0.0,0.0,-dom.GetParticle(-4)->Props.m*9.8;
    dom.GetParticle(-4)->FixVeloc();
    dom.GetParticle(-4)->vzf=false;

    
    dom.Solve(20.0,dt,0.1,NULL,NULL,"drillb",true);
    return 0;
}
MECHSYS_CATCH
