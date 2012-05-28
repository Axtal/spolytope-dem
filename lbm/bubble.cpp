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

using std::cout;
using std::endl;
struct UserData
{
    double            Orig;
    double             Amp;
    double             ome;
    double              Tf;
    Vec3_t               g;
    Array<Cell *>   Center;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double rho;
    //rho = 1.0 + 0.001*0.5*(1-cos(dat.ome*dom.Time));
    //rho = 1.0 + 0.1*sin(dat.ome*dom.Time);
    rho = dat.Orig + dat.Amp*sin(dat.ome*dom.Time);
    for(size_t i=0;i<dat.Center.Size();i++)
    {
        dat.Center[i]->Initialize(rho,OrthoSys::O);
    }
}


int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double Gmix  = 0.001;
    double Gs0   = 0.001;
    double Gs1   = 0.001;
    double R     = 10.0;
    double rho0  = 1.0;
    double rho1  = 0.01;
    double ome   = 1.0;
    double Orig  = 1.0;
    double Amp   = 1.0;
    double Tf    = 5000.0;

    if (argc>=2)
    {
        Nproc  =atoi(argv[ 1]);
        Gmix   =atof(argv[ 2]);
        Gs0    =atof(argv[ 3]);
        Gs1    =atof(argv[ 4]);
        R      =atof(argv[ 5]);
        rho0   =atof(argv[ 6]);
        rho1   =atof(argv[ 7]);
        ome    =atof(argv[ 8]);
        Orig   =atof(argv[ 9]);
        Amp    =atof(argv[10]);
        Tf     =atof(argv[11]);
    }


    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    size_t nx   = 200, ny = 100;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    //dat.g           = 0.0001,-0.001,0.0;
    //dat.g           = 0.0,-0.001,0.0;
    //dat.g          *= 0.5*ny/100.0;
    //dat.g           = 0.0,-0.0005,0.0;
    dat.ome         = 2*M_PI*ome/Tf;
    dat.Tf          = Tf;
    dat.Amp         = Amp;
    dat.Orig        = Orig;



    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    //for (size_t i=0;i<ny;i++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    //}

    // Set inner drop
    int obsX = nx/2, obsY = ny/2;
    int r1 =  0.1*R;
    int r2 =  R;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(r1,2.0)) // circle equation
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(rho0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(rho1,V);
            dat.Center.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0)));
		}
		else if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(r2,2.0)) 
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(rho0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(rho1,V);
		}
        else
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(rho1,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(rho0,V);
        }
    }

    // Set parameters
    Dom.Lat[0].G =  0.0 ;
    Dom.Lat[0].Gs=  Gs0;
    Dom.Lat[1].G =  0.0 ;
    Dom.Lat[1].Gs=  Gs1  ;
    Dom.Gmix     =  Gmix ;
    //Dom.Lat[0].G = -150.0;
    //Dom.Lat[0].Gs=  400.0;
    //Dom.Lat[1].G = -200.0;
    //Dom.Lat[1].Gs= -400.0;
    //Dom.Gmix     =  0.0001;
    //
    //Dom.Lat[0].G =  0.0;
    //Dom.Lat[0].Gs=  2.0;
    //Dom.Lat[1].G =  0.0;
    //Dom.Lat[1].Gs= -2.0;
    //Dom.Gmix     =  4.0;

    Dom.Solve(Tf,0.01*Tf,Setup,NULL,"bubble",true,Nproc);


    return 0;
}
MECHSYS_CATCH

