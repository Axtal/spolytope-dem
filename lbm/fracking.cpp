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
    double              Tf;
    Array<Cell *>   Center;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double rho;
    rho = dat.Orig + dom.Time*(dat.Amp - dat.Orig)/dat.Tf;
    for(size_t i=0;i<dat.Center.Size();i++)
    {
        dat.Center[i]->Initialize(rho,OrthoSys::O);
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double Amax  = 10.0;
    double nu    = 0.2;
    double Orig  = 1.0;
    double Amp   = 1.0;
    double Tf    = 5000.0;
    double dt    = 1.0;
    double Alpha = 10.0;
    double r1    = 2.5;
    double r2    = 5.0;
    double Kn    = 300.0;
    double Kt    = 300.0;
    double Gn    = 1.0;
    double Gt    = 0.0;
    double Mu    = 0.3;
    double Bn    = 10.0;
    double Bt    = 10.0;
    double Bm    = 0.0;
    double Eps   = 0.01;

    if (argc>=2)
    {
        Nproc  =atoi(argv[ 1]);
        Amax   =atof(argv[ 2]);
        Orig   =atof(argv[ 3]);
        Amp    =atof(argv[ 4]);
        Tf     =atof(argv[ 5]);
        dt     =atof(argv[ 6]);
        Alpha  =atof(argv[ 6]);
    }


    size_t nx   = 500, ny = 500, nz = 3;

    r2 = nx/20;
    r1 = 0.8*r2;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), 1.0, dt);
    UserData dat;
    Dom.UserData = &dat;
    dat.Tf          = Tf;
    dat.Amp         = Amp;
    dat.Orig        = Orig;
    Dom.Alpha       = Alpha;

    size_t n_divisions = 10;

    Mesh::Unstructured mesh(/*NDim*/2);
    mesh.Set    (4+n_divisions,4+n_divisions, 1, 1);                            // 8 points, 8 segments, 1 region, 1 hole
    mesh.SetReg (0, -1, Amax, 0.2, 0.2);                 // id, tag, max{area}, x, y <<<<<<< regions
    mesh.SetHol (0, 0.5*nx, 0.5*ny);                     // id, x, y <<<<<<< holes
    mesh.SetPnt (0, -1, 0.0, 0.0);                       // id, vtag, x, y <<<<<< points
    mesh.SetPnt (1, -2,  nx, 0.0);                       // id, vtag, x, y
    mesh.SetPnt (2, -3,  nx,  ny);                       // id, vtag, x, y
    mesh.SetPnt (3, -4, 0.0,  ny);                       // id, vtag, x, y
    mesh.SetSeg (0, -10,  0, 1);                         // id, etag, L, R <<<<<<<<<<<< segments
    mesh.SetSeg (1, -20,  1, 2);                         // id, etag, L, R
    mesh.SetSeg (2, -30,  2, 3);                         // id, etag, L, R
    mesh.SetSeg (3, -40,  3, 0);                         // id, etag, L, R
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetPnt(i+4, 0, r2*cos(2*i*M_PI/n_divisions)+0.5*nx,r2*sin(2*i*M_PI/n_divisions)+0.5*ny);
    }
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetSeg(i+4, 0, i + 4, (i+1)%n_divisions + 4);
    }

    mesh.Generate ();
    Dom.GenFromMesh(mesh,/*spheroradius*/1.0,/*density*/3.0,/*iscohesive*/false,/*montecarlo mass properties*/false,/*thickness*/5.0);
    Dom.Center(Vec3_t(0.5*(nx-1),0.5*(ny-1),0.5*(nz-1)));
	
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt Eps Kn Kt Mu"   ,Bn,Bt,Bm,Gn,Gt,Eps,Kn,Kt,Mu);
    Dom.SetProps(B);
    
    for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
	for (size_t k=0; k<nz; ++k)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,OrthoSys::O);
		if (pow((int)(i)-0.5*nx,2.0) + pow((int)(j)-0.5*ny,2.0) <= pow(r1,2.0)) // circle equation
		{
            dat.Center.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,k)));
		}
    }

    Dom.Solve(Tf,0.01*Tf,Setup,NULL,"fracking",true,Nproc);


    return 0;
}
MECHSYS_CATCH

