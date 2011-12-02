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
    double Tf    = 4.5;
    double dt    = 1.0e-4;
    double dtOut = 0.01*Tf;
    double ome   = 0.0*2*M_PI/Tf;

    double Kn = 1.0e5;        // Normal stiffness
    double Kt = 1.0e5;        // Tangential stiffness
    double Gn = 8.0;          // Normal dissipative coefficient
    double Gt = 8.0;          // Tangential dissipative coefficient
    double Mu = 0.3;          // Microscopic friction coefficient
    double Bn = 2.0e5;        // Cohesion normal stiffness
    double Bt = 2.0e5;        // Cohesion tangential stiffness
    double Bm = 2.0e2;        // Cohesion torque stiffness
    double Eps = 0.05;        // Threshold for breking bonds

    // domain and User data
    DEM::Domain dom;
    dom.Alpha  = 0.05;
    dom.CamPos = Vec3_t(0.0,20.0,0.0);

    // particle
    
    size_t Ind = dom.Particles.Size();
    Ind = dom.Particles.Size();
    dom.AddVoroPack (-1, 0.1, 4,4,4, 4,4,4, 3.0, true, true, 1000, 1.0);
    for (size_t i=Ind;i<dom.Particles.Size();i++)
    {
        Vec3_t trans(0.0,0.0,5.0);
        dom.Particles[i]->Translate(trans);
    }
    Ind = dom.Particles.Size();
    dom.AddVoroPack (-1, 0.1, 4,4,4, 4,4,4, 3.0, true, true, 3000, 1.0);
    for (size_t i=Ind;i<dom.Particles.Size();i++)
    {
        Vec3_t trans(0.0,0.0,10.0);
        dom.Particles[i]->Translate(trans);
    }
    Ind = dom.Particles.Size();
    dom.AddVoroPack (-1, 0.1, 4,4,4, 4,4,4, 3.0, true, true, 5000, 1.0);
    for (size_t i=Ind;i<dom.Particles.Size();i++)
    {
        Vec3_t trans(0.0,0.0,15.0);
        dom.Particles[i]->Translate(trans);
    }
    
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Initialize(i);
        dom.Particles[i]->Ff = 0.0,0.0,-dom.Particles[i]->Props.m*9.8;
    }

    dom.AddRice(-2,Vec3_t( 10.0,0.0,0.0),8.0,10.0,3.0,M_PI/2.0,&OrthoSys::e0);
    dom.Particles[dom.Particles.Size()-1]->FixVeloc();
    dom.Particles[dom.Particles.Size()-1]->w = Vec3_t(0.0,0.0, ome);
    dom.AddRice(-3,Vec3_t(-10.0,0.0,0.0),8.0,10.0,3.0,M_PI/2.0,&OrthoSys::e0);
    dom.Particles[dom.Particles.Size()-1]->FixVeloc();
    dom.Particles[dom.Particles.Size()-1]->w = Vec3_t(0.0,0.0,-ome);
    dom.AddRice(-4,Vec3_t( 10.0,0.0,0.0),1.0, 5.0,3.0,M_PI/2.0,&OrthoSys::e1);
    dom.Particles[dom.Particles.Size()-1]->FixVeloc();
    dom.Particles[dom.Particles.Size()-1]->w = Vec3_t(0.0,-ome, 0.0);
    dom.AddRice(-5,Vec3_t(-10.0,0.0,0.0),1.0, 5.0,3.0,M_PI/2.0,&OrthoSys::e1);
    dom.Particles[dom.Particles.Size()-1]->FixVeloc();
    dom.Particles[dom.Particles.Size()-1]->w = Vec3_t(0.0, ome, 0.0);

    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu,Bn,Bt ,Bm ,     Eps);
    B.Set(-2,"Kn Kt Gn Gt Mu Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-3,"Kn Kt Gn Gt Mu Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B);

    dom.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, "grinding", true);

    return 0;
}
MECHSYS_CATCH

