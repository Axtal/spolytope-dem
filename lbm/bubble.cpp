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
    double             ome;
    double              Tf;
    Vec3_t               g;
    Array<Cell *>   Center;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t j=0;j<dom.Lat.Size();j++)
    for (size_t i=0;i<dom.Lat[j].Cells.Size();i++)
    {
        Cell * c = dom.Lat[j].Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
    double rho;
    rho = 800.0 + 300.0*0.5*(1-cos(dat.ome*dom.Time));
    //rho = 1600.0;
    //rho = 1.0 + 0.1*0.5*(1-cos(dat.ome*dom.Time));
    //if (dom.Time<0.5*dat.Tf) rho = 1370.0;
    //else                     rho =  900.0;     
    for(size_t i=0;i<dat.Center.Size();i++)
    {
        dat.Center[i]->Initialize(rho,OrthoSys::O);
    }
}


int main(int argc, char **argv) try
{
    //double Gs = atof(argv[1]); 
    double Gs = 300;


    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    size_t nx   = 200, ny = 100;
    double Tf   = 10000.0;
    double ome  = 1.0;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    //dat.g           = 0.0001,-0.001,0.0;
    //dat.g           = 0.0,-0.001,0.0;
    //dat.g          *= 0.5*ny/100.0;
    dat.g           = 0.0,-0.0005,0.0;
    dat.ome         = 2*M_PI*ome/Tf;
    dat.Tf          = Tf;



    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }

    // Set inner drop
    int obsX = nx/2, obsY = ny/2, obsZ = 3*ny/4;
    //int radius =  ny/4.0;
    int r1 =  ny/20.0;
    int r2 =  ny/3.0;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(r1,2.0)) // circle equation
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(800.0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1.0 ,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.001,V);
            dat.Center.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0)));
		}
		else if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsZ,2.0) <= pow(r2,2.0)) 
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(800.0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1  ,V);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.001,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(1.0  ,V);
		}
        else
        {
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,V);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.001,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(1.0  ,V);
        }
    }

    // Set parameters
    Dom.Lat[0].G = -150.0;
    Dom.Lat[0].Gs= -300.0;
    Dom.Lat[1].G = -200.0;
    Dom.Lat[1].Gs= -Gs;
    Dom.Gmix     =  0.0001;
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

    Dom.Solve(Tf,0.01*Tf,Setup,NULL,"bubble");


    return 0;
}
MECHSYS_CATCH

