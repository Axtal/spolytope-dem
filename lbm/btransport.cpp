/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <iostream>
#include <stdlib.h>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;


struct UserData
{
    std::ofstream oss_st;       ///< file for surface tension calculation
    double          rho0;
    double          rho1;
    Vec3_t             g;
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.g;
    }
}

void Report (LBM::Domain & dom, void * UD)
{
}



int main(int argc, char **argv) try
{
    size_t Nproc = 1;
    size_t nx    = 200;
    size_t ny    = 300;
    double g     = 0.1;
    double Gmix  = 0.001;
    double Gs0   = -200.0;
    double Gs1   = 0.0;
    double R     = 10.0;
    double rho0  = 1.0;
    double rho1  = 0.01;
    double Tf    = 5000.0;
    double dtout = 50.0;
    double por   = 0.4;
    if (argc>=2)
    {
        Nproc  =atoi(argv[1]);
        nx     =atoi(argv[2]);
        ny     =atoi(argv[3]);
        g      =atof(argv[4]);
        Gmix   =atof(argv[5]);
        Gs0    =atof(argv[6]);
        Gs1    =atof(argv[7]);
        R      =atof(argv[8]);
        rho0   =atof(argv[9]);
        rho1   =atof(argv[10]);
        Tf     =atof(argv[11]);
        dtout  =atof(argv[12]);
        por    =atof(argv[13]);
    }
    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.rho0     = rho0;
    dat.rho1     = rho1;
    dat.g        = Vec3_t(0.0,g,0.0);

    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->Gs      = 0.0;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->Gs      = 0.0;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->Gs      = 0.0;
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->Gs      = 0.0;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->Gs      = 0.0;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->Gs      = 0.0;
        Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->Gs      = 0.0;
        Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->Gs      = 0.0;
    }

    srand(1000);
    size_t ntries = 0;
    double n      = 1.0/6.0;
    Array<double> Radii;
    Array<Vec3_t> Xs;
    double        rc = 4.0;
    ntries = 0;
    
    n      = 5.0/6.0;
    Radii.Resize(0);
    Xs.Resize(0);

    while (1-Dom.Lat[0].SolidFraction()/n>por)
    {
        ntries++;
        if (ntries>1.0e4) throw new Fatal("Too many tries to achieved requested porosity, please increase it");
        double Rmax = 0.1;
        double Rmin = 0.4*Rmax;
        double r  = ((Rmin*Rmax/(Rmax - double(rand())/RAND_MAX*(Rmax - Rmin))))*nx;
        double DY = 1.0/6.0*ny;
        double yc = DY + (ny - DY)*double(rand())/RAND_MAX;
        double xc = nx*double(rand())/RAND_MAX;
        Vec3_t X(xc,yc,0.0);
        bool invalid = false;
        for (size_t i=0;i<Radii.Size();i++)
        {
            if (norm(Xs[i] - X) < r + Radii[i] + rc)
            {
                invalid = true;
                break;
            }
        }
        if (invalid) continue;
        Dom.Lat[0].SolidDisk(Vec3_t(xc,yc,0.0),r);
        Dom.Lat[1].SolidDisk(Vec3_t(xc,yc,0.0),r);
        Radii.Push(r);
        Xs.Push(Vec3_t(xc,yc,0.0));
    }

    int obsX = nx/2, obsY = ny/12;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.001*rho1,OrthoSys::O);
        Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.999*rho1,OrthoSys::O);
        for (obsX = 1.5*R;obsX<nx-1.5*R;obsX+=3.0*R)
        {
		    if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(R,2.0)) // circle equation
		    {
                Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.999*rho0,OrthoSys::O);
                Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.001*rho0,OrthoSys::O);
		    }
        }
    }

    // Set parameters
    Dom.Gmix      = Gmix;
    Dom.Lat[0].Gs = Gs0;
    Dom.Lat[1].Gs = Gs1;

    Dom.Solve(Tf,dtout,Setup,Report,"btransport",true,Nproc);

    return 0;
}
MECHSYS_CATCH
