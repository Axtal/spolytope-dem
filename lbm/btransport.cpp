/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2013 Sergio Torres                                     *
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
    std::ofstream  oss_st;       ///< file for surface tension calculation
    double            rho;
    double            inj;
    Vec3_t              g;
    Array<Cell *>      c0;
    Array<Cell *>      c1;
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat[0].Ncells;i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.g;
    }
    for (size_t i=0;i<dat.c0.Size();i++)
    {
        //dat.c0[i]->Initialize(0.999*(0.99 + (0.02*rand())/RAND_MAX)*dat.inj,OrthoSys::O);
        //dat.c1[i]->Initialize(0.001*(0.99 + (0.02*rand())/RAND_MAX)*dat.inj,OrthoSys::O);
        dat.c0[i]->Initialize(0.999*dat.inj,OrthoSys::O);
        dat.c1[i]->Initialize(0.001*dat.inj,OrthoSys::O);
    }
}

void Report (LBM::Domain & dom, void * UD)
{
}

int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc>=3) Nproc = atoi(argv[2]);

    size_t nx    = 200;
    size_t ny    = 300;
    double g     = 0.1;
    double Gmix  = 0.001;
    double Gs0   = -200.0;
    double Gs1   = 0.0;
    double R     = 10.0;
    double space = 10.0;
    size_t Ninj  = 4;
    double Rmin  = 10.0;
    double Rmax  = 10.0;
    double rho  = 1.0;
    double inj  = 1.1;
    double Tf    = 5000.0;
    double dtOut = 50.0;
    double por   = 0.5;
    size_t seed  = 1000;


    infile >> filename;  infile.ignore(200,'\n');
    infile >> nx;        infile.ignore(200,'\n');
    infile >> ny;        infile.ignore(200,'\n');
    infile >> g;         infile.ignore(200,'\n');
    infile >> Gmix;      infile.ignore(200,'\n');
    infile >> Gs0;       infile.ignore(200,'\n');
    infile >> Gs1;       infile.ignore(200,'\n');
    infile >> R;         infile.ignore(200,'\n');
    infile >> space;     infile.ignore(200,'\n');
    infile >> Ninj;      infile.ignore(200,'\n');
    infile >> Rmin;      infile.ignore(200,'\n');
    infile >> Rmax;      infile.ignore(200,'\n');
    infile >> rho;       infile.ignore(200,'\n');
    infile >> inj;       infile.ignore(200,'\n');
    infile >> Tf;        infile.ignore(200,'\n');
    infile >> dtOut;     infile.ignore(200,'\n');
    infile >> por;       infile.ignore(200,'\n');
    infile >> seed;      infile.ignore(200,'\n');

    Array<double> XC;
    Array<double> YC;
    Array<double> RC;
    double xmax = 0.0;
    double ymax = 0.0;

    bool fileexist = true;
    if (!Util::FileExists(filename))
    {
        fileexist = false;
    }
    else 
    {
        std::ifstream infile2(filename.CStr());
        while (!infile2.eof())
        {
            double xc,yc,rc;
            infile2 >> xc;
            infile2 >> yc;
            infile2 >> rc;
            XC.Push(xc);
            YC.Push(yc);
            RC.Push(rc);
            if (xc+rc>xmax) xmax = xc+rc;
            if (yc+rc>ymax) ymax = yc+rc;
        }
        ny = Rmin*ymax;
        nx = Rmin*xmax;
    }

    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.rho      = rho;
    dat.inj      = inj;
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

    srand(seed);
    size_t ntries = 0;
    double n      = 1.0/6.0;
    Array<double> Radii;
    Array<Vec3_t> Xs;
    double        rc = 4.0;
    ntries = 0;
    
    n      = 5.0/6.0;
    Radii.Resize(0);
    Xs.Resize(0);

    if (fileexist)
    {
        for (size_t i=0;i<XC.Size();i++)
        {
            if (Rmin*(YC[i]-RC[i])>Rmax)
            {
                Dom.Lat[0].SolidDisk(Vec3_t(Rmin*XC[i],Rmin*YC[i],0.0),Rmin*por*RC[i]);
                Dom.Lat[1].SolidDisk(Vec3_t(Rmin*XC[i],Rmin*YC[i],0.0),Rmin*por*RC[i]);
            }
        }
    }
    else
    {
        while (1-Dom.Lat[0].SolidFraction()/n>por)
        {
            ntries++;
            if (ntries>1.0e4) throw new Fatal("Too many tries to achieved requested porosity, please increase it");
            //double Rmax = 0.1;
            //double Rmin = 0.4*Rmax;
            double r  = ((Rmin*Rmax/(Rmax - double(rand())/RAND_MAX*(Rmax - Rmin))));
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
    }

    int obsX = nx/2, obsY = (fileexist) ? Rmax/2 : ny/12;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.001*rho,OrthoSys::O);
        Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.999*rho,OrthoSys::O);
        for (int obsx = obsX-0.5*(Ninj-1)*space;obsx<=obsX+0.5*(Ninj-1)*space;obsx+=space)
        {
		    if (pow((int)(i)-obsx,2.0) + pow((int)(j)-obsY,2.0) <= pow(R,2.0)) // circle equation
		    {
                Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.999*inj,OrthoSys::O);
                Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.001*inj,OrthoSys::O);
                dat.c0.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0)));
                dat.c1.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,0)));
		    }
        }
    }

    // Set parameters
    Dom.Gmix      = Gmix;
    Dom.Lat[0].Gs = Gs0;
    Dom.Lat[1].Gs = Gs1;

    //Dom.WriteXDMF("btransport");

    Dom.Solve(Tf,dtOut,Setup,Report,"btransport",true,Nproc);

    return 0;
}
MECHSYS_CATCH

