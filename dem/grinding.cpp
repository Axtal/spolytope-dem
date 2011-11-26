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
    // domain and User data
    DEM::Domain dom;
    dom.Alpha  = 0.05;
    dom.CamPos = Vec3_t(0.1*Lx, 0.7*(Lx+Ly+Lz), 0.15*Lz); // position of camera

    // particle
    {
        size_t Ind = dom.Particles.Size();
        dom.AddVoroPack (-1, R, 4,4,4, 4,4,4, rho, Cohesion, Cohesion, seed, fraction);
        for (size_t i=Ind;i<dom.Particles.Size();i++)
        {
            Vec3_t trans(-5.0,0.0,0.0);
            dom.Particles[i]->Translate(trans);
        }
        Ind = dom.Particles.Size();
        dom.AddVoroPack (-1, R, 4,4,4, 4,4,4, rho, Cohesion, Cohesion, seed, fraction);
        for (size_t i=Ind;i<dom.Particles.Size();i++)
        {
            Vec3_t trans( 5.0,0.0,0.0);
            dom.Particles[i]->Translate(trans);
        }
        Ind = dom.Particles.Size();
        dom.AddVoroPack (-1, R, 4,4,4, 4,4,4, rho, Cohesion, Cohesion, seed, fraction);
        for (size_t i=Ind;i<dom.Particles.Size();i++)
        {
            Vec3_t trans(0.0,0.0,5.0);
            dom.Particles[i]->Translate(trans);
        }
        Ind = dom.Particles.Size();
        dom.AddVoroPack (-1, R, 4,4,4, 4,4,4, rho, Cohesion, Cohesion, seed, fraction);
        for (size_t i=Ind;i<dom.Particles.Size();i++)
        {
            Vec3_t trans(0.0,0.0,-5.0);
            dom.Particles[i]->Translate(trans);
        }
    }
    

    // properties of particles prior the triaxial test
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    B.Set(-2,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-3,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-6,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B.Set(-7,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt,Mu ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B);

    // stage 1: isotropic compresssion  //////////////////////////////////////////////////////////////////////
    String fkey_a(filekey+"_a");
    String fkey_b(filekey+"_b");
    Vec3_t  sigf;                      // final stress state
    bVec3_t peps(false, false, false); // prescribed strain rates ?
    Vec3_t  depsdt(0.0,0.0,0.0);       // strain rate

    sigf =  Vec3_t(-p0,-p0,-p0);
    ResetEps  (dom,dat);
    SetTxTest (sigf, peps, depsdt,0,0,false,dat,dom);
    dat.tspan = T0/2.0 - dom.Time;
    dom.Solve  (/*tf*/T0/2.0, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_a.CStr(),RenderVideo);
    SetTxTest (sigf, peps, depsdt,0,0,false,dat,dom);
    dat.tspan = T0 - dom.Time;
    dom.Solve (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_b.CStr(),RenderVideo);

    // stage 2: The proper triaxial test /////////////////////////////////////////////////////////////////////////
    String fkey_c(filekey+"_c");
    Vec3_t lf;
    pqth2L (pf, qf, thf, lf, "cam");
    sigf   = lf(0), lf(1), lf(2);
    peps   = bVec3_t(pssrx, pssry, pssrz);
    depsdt = Vec3_t(srx/(Tf-dom.Time), sry/(Tf-dom.Time), srz/(Tf-dom.Time));
    
    // run
    ResetEps  (dom,dat);
    SetTxTest (sigf, peps, depsdt, thf*M_PI/180, alpf*M_PI/180, isfailure, dat, dom);
    dat.tspan = Tf - dom.Time;
    dom.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup, &Report, fkey_c.CStr(),RenderVideo);

    return 0;
}
MECHSYS_CATCH

