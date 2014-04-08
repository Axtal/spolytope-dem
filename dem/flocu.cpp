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

struct UserData
{
    Array<Particle*> Left;
    Array<Particle*> Right;
};

void Setup1 (DEM::Domain & dom, void * UD)
{
    double G = 0.001;

    for (size_t i = 0; i < dom.Particles.Size()-1; i++)
    {
        for (size_t j = i+1; j < dom.Particles.Size() ; j++)
        {
            if(dom.Particles[i]->IsFree()&&dom.Particles[j]->IsFree())
            {
                Vec3_t x1 = dom.Particles[i]->x;
                Vec3_t x2 = dom.Particles[j]->x;
                Vec3_t F  = (G/(pow(norm(x2-x1),3)))*(x2-x1);
                dom.Particles[i]->F += F;
                dom.Particles[j]->F -= F;
            }
        }
    }
}

void Setup2 (DEM::Domain & dom, void * UD)
{
    double G = 0.001;
    double K = 0.001;
    Vec3_t Xcmin,Xcmax;
    dom.BoundingBox(Xcmin, Xcmax);
    Vec3_t Xc = 0.5*(Xcmin + Xcmax);
    for (size_t i = 0; i < dom.Particles.Size()-1; i++)
    {
        for (size_t j = i+1; j < dom.Particles.Size() ; j++)
        {
            if(dom.Particles[i]->IsFree()&&dom.Particles[j]->IsFree())
            {
                Vec3_t x1 = dom.Particles[i]->x;
                Vec3_t x2 = dom.Particles[j]->x;
                Vec3_t F  = (G/(pow(norm(x2-x1),3)))*(x2-x1);
                dom.Particles[i]->F += F;
                dom.Particles[j]->F -= F;
            }
        }
        dom.Particles[i]->F+=K*((dom.Particles[i]->x(0) - Xc(0)<0)?-1:1)*OrthoSys::e0;
    }
}

void Setup3 (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double G = 0.001;
    for (size_t i = 0; i < dat.Left.Size()-1; i++)
    {
        for (size_t j = i+1; j < dat.Left.Size() ; j++)
        {
            if(dat.Left[i]->IsFree()&&dat.Left[j]->IsFree())
            {
                Vec3_t x1 = dat.Left[i]->x;
                Vec3_t x2 = dat.Left[j]->x;
                Vec3_t F  = (G/(pow(norm(x2-x1),3)))*(x2-x1);
                dat.Left[i]->F += F;
                dat.Left[j]->F -= F;
            }
        }
    }
    for (size_t i = 0; i < dat.Right.Size()-1; i++)
    {
        for (size_t j = i+1; j < dat.Right.Size() ; j++)
        {
            if(dat.Right[i]->IsFree()&&dat.Right[j]->IsFree())
            {
                Vec3_t x1 = dat.Right[i]->x;
                Vec3_t x2 = dat.Right[j]->x;
                Vec3_t F  = (G/(pow(norm(x2-x1),3)))*(x2-x1);
                dat.Right[i]->F += F;
                dat.Right[j]->F -= F;
            }
        }
    }
}

void Report (DEM::Domain & dom, void *UD)
{
}

int main(int argc, char **argv) try
{
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = 0.01;

    // particle
    double Lx       = 10.0;
    size_t nx       = 10;
    double rho      = 3.0;
    size_t seed     = 1000;
    double fraction = 0.1;
    double Tf       = 300.0;
    double dt       = 1.0e-3;
    double dtOut    = 1.0;
    double v0       = 1.0;
    

    dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction);
    for (size_t i = 0; i < dom.Particles.Size(); i++)
    {
        dom.Particles[i]-> v = v0*Vec3_t((1.0*rand())/RAND_MAX-0.5,(1.0*rand())/RAND_MAX-0.5,(1.0*rand())/RAND_MAX-0.5);
    }
    Dict B;
    B.Set(-1,"Mu ",0.0);
    dom.SetProps(B);
    dom.CamPos = Vec3_t(0.1*Lx, 0.7*(3.0*Lx), 0.15*Lx); // position of camera
    dom.GenBoundingBox (/*InitialTag*/-2, 0.2, /*Cf*/1.0);
    dom.GetParticle (-2)->FixVeloc();
    dom.GetParticle (-3)->FixVeloc();
    dom.GetParticle (-4)->FixVeloc();
    dom.GetParticle (-5)->FixVeloc();
    dom.GetParticle (-6)->FixVeloc();
    dom.GetParticle (-7)->FixVeloc();
        
    dom.Solve(/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup1, &Report, "flocu1", true);
    Vec3_t Xcmin,Xcmax;
    dom.BoundingBox(Xcmin, Xcmax);
    Vec3_t Xc = 0.5*(Xcmin + Xcmax);
    for (size_t i = 0; i < dom.Particles.Size(); i++)
    {
        if (dom.Particles[i]->IsFree())
        {
            if (dom.Particles[i]->x(0)<Xc(0)) dat.Left.Push(dom.Particles[i]);
            else                              dat.Right.Push(dom.Particles[i]);
        }
    }
    dom.Solve(/*tf*/2*Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup3, &Report, "flocu2", true);

    return 0;
}
MECHSYS_CATCH
