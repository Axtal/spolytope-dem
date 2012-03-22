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
// Vorticity test

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>

struct UserData
{
    double          Kn;
    Vec3_t           g;
    Vec3_t        Xmin;
    Vec3_t        Xmax;
    double          Tf;
    double          dt;
    double         DPz;
};

void Setup (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        //dom.Particles[i]->Ff = dom.Particles[i]->Props.m*dat.g - 0.8*dom.Particles[i]->Props.m*dom.Particles[i]->v;
        dom.Particles[i]->Ff = dom.Particles[i]->Props.m*dat.g;
        double delta;
        delta =   dat.Xmin(0) - dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(0) += dat.Kn*delta;
        delta = - dat.Xmax(0) + dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(0) -= dat.Kn*delta;
        delta =   dat.Xmin(1) - dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) + dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) -= dat.Kn*delta;
        delta =   dat.Xmin(2) - dom.Particles[i]->x(2) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(2) += dat.Kn*delta;
        //delta = - dat.Xmax(2) + dom.Particles[i]->X(2) + dom.Particles[i]->R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(2) -= dat.Kn*delta;
    }
}

void Report (DEM::Domain & dom, void * UD)
{
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc==3) Nproc = atoi(argv[2]);

    bool   Render   = true;
    size_t seed     = 100;
    double fraction = 0.8;
    double Rmin     = 0.9;
    double nu       = 0.1;
    size_t nx       = 100;
    size_t ny       = 100;
    size_t nz       = 100;
    double g        = 0.001;
    double Tf       = 4000.0;
    double dtOut    = 20.0;
    double dx       = 1.0;
    double dt       = 1.0;
    double Kn       = 1.0e3;
    double DPz      = 0.0;
    double R1       = 2.0;
    double R2       = 20.0;

    infile >> Render;    infile.ignore(200,'\n');
    infile >> seed;      infile.ignore(200,'\n');
    infile >> fraction;  infile.ignore(200,'\n');
    infile >> Rmin;      infile.ignore(200,'\n');
    infile >> nu;        infile.ignore(200,'\n');
    infile >> nx;        infile.ignore(200,'\n');
    infile >> ny;        infile.ignore(200,'\n');
    infile >> nz;        infile.ignore(200,'\n');
    infile >> g;         infile.ignore(200,'\n');
    infile >> Tf;        infile.ignore(200,'\n');
    infile >> dtOut;     infile.ignore(200,'\n');
    infile >> dx;        infile.ignore(200,'\n');
    infile >> dt;        infile.ignore(200,'\n');
    infile >> Kn;        infile.ignore(200,'\n');
    infile >> DPz;       infile.ignore(200,'\n');
    infile >> R1;        infile.ignore(200,'\n');
    infile >> R2;        infile.ignore(200,'\n');

    DEM::Domain Dom;
    UserData dat;
    Dom.UserData = &dat;
    dat.g        = 0.0,0.0,-g;
    dat.Xmin     = 0.0,0.0,dx*nz/10.0;
    dat.Xmax     = nx*dx,ny*dx,nz*dx;
    dat.Kn       = Kn;
    dat.Tf       = Tf;
    dat.dt       = dt;
    dat.DPz      = DPz;

    Dom.Alpha    = std::min(R1,R2);
    

    Dom.GenSpheresBox (-1, dat.Xmin , dat.Xmax - Vec3_t(0.0,0.0,0.5*nz*dx), /*R*/R1, 2.5, seed, fraction, Rmin); ///< Create an array of spheres
    Dom.GenSpheresBox (-2, dat.Xmax - Vec3_t(nx*dx,ny*dx,0.5*nz*dx), dat.Xmax + Vec3_t(0.0,0.0,0.3*nz*dx), /*R*/R2, 2.5, seed, fraction, Rmin); ///< Create an array of spheres
    Dom.CamPos = Vec3_t(nz*dx,nz*dx,nx*dx);
    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        Dom.Particles[i]->Props.Kn = Kn;
        Dom.Particles[i]->Props.Kt = Kn;
    }

    Dom.Solve(Tf,dt,dtOut,&Setup,&Report,filekey.CStr(),Render,Nproc);
    Dom.Save(filekey.CStr());
    return 0;
}
MECHSYS_CATCH

