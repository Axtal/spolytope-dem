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
    bool               rev;
    double            Orig;
    double         Gsolid0;
    double         Gsolid1;
    double             Amp;
    double             ome;
    double              Tf;
    double          revrad1;
    double          revrad2;
    Vec3_t               g;
    Array<Cell *>   Center;
    Array<size_t>      Top;
    Array<size_t>   Bottom;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
    {
        dom.Lat[0].Cells[i]->BForcef = dom.Lat[0].Cells[i]->Rho*dat.g;
    }
    double rho;
    //rho = 1.0 + 0.001*0.5*(1-cos(dat.ome*dom.Time));
    //rho = 1.0 + 0.1*sin(dat.ome*dom.Time);
    rho = dat.Orig + dat.Amp*sin(dat.ome*dom.Time);
    for(size_t i=0;i<dat.Center.Size();i++)
    {
        dat.Center[i]->Initialize(rho,OrthoSys::O);
    }
    if (dat.rev)
    {
        for (size_t i=0;i<dat.Top.Size();i++)
        {
            Cell * c0  = dom.Lat[0].Cells[dat.Top[i]];
            Cell * c1  = dom.Lat[1].Cells[dat.Top[i]];
            Cell * nb  = dom.Lat[0].Cells[c0->Neighs[4]];
            Cell * nbr = dom.Lat[0].Cells[c0->Neighs[7]];
            Cell * nbl = dom.Lat[0].Cells[c0->Neighs[8]];
            nbr = dom.Lat[0].GetCell(iVec3_t(std::min(size_t(c0->Index(0) + dat.revrad1),size_t(dom.Lat[0].Ndim(0))),c0->Index(1)-1,c0->Index(2)));
            nbl = dom.Lat[0].GetCell(iVec3_t(std::max(size_t(c0->Index(0) - dat.revrad1),size_t(0)                 ),c0->Index(1)-1,c0->Index(2)));
            if (nb->Rho>1.0&&(nbr->Rho>1.0&&nbl->Rho>1.0)&&c0->Gs>0.0) 
            {
                //std::cout << dom.Time << " " << c0->Index << std::endl;
                c0->Gs = -1.0;
                c1->Gs = -1.0;
            }
            else
            { 
                nbr = dom.Lat[0].GetCell(iVec3_t(std::min(size_t(c0->Index(0) + dat.revrad2),size_t(dom.Lat[0].Ndim(0))),c0->Index(1)-1,c0->Index(2)));
                nbl = dom.Lat[0].GetCell(iVec3_t(std::max(size_t(c0->Index(0) - dat.revrad2),size_t(0)                 ),c0->Index(1)-1,c0->Index(2)));
                if (c0->Gs<0.0&&(nbr->Rho<1.0&&nbl->Rho<1.0)&&nb->Rho<1.0)
                {
                    //std::cout << dom.Time << " " << c0->Index << std::endl;
                    c0->Gs =  1.0;
                    c1->Gs =  1.0;
                }
            }
        }


        //dom.Lat[0].Gs = dat.Gsolid0*cos(dat.ome*dom.Time);
        //dom.Lat[1].Gs = dat.Gsolid1*cos(dat.ome*dom.Time);
    }
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc==3) Nproc = atoi(argv[2]);

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
    bool   rev   = false;
    double revrad1= 4.0;
    double revrad2= 4.0;
    
    infile >>    Gmix;        infile.ignore(200,'\n');
    infile >>    Gs0;         infile.ignore(200,'\n');
    infile >>    Gs1;         infile.ignore(200,'\n');
    infile >>    R;           infile.ignore(200,'\n');
    infile >>    rho0;        infile.ignore(200,'\n');
    infile >>    rho1;        infile.ignore(200,'\n');
    infile >>    ome;         infile.ignore(200,'\n');
    infile >>    Orig;        infile.ignore(200,'\n');
    infile >>    Amp;         infile.ignore(200,'\n');
    infile >>    Tf;          infile.ignore(200,'\n');
    infile >>    rev;         infile.ignore(200,'\n');
    infile >>    revrad1;     infile.ignore(200,'\n');
    infile >>    revrad2;     infile.ignore(200,'\n');
    

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
    dat.g           = 0.0,0.0001,0.0;
    dat.ome         = 2*M_PI*ome/Tf;
    dat.Tf          = Tf;
    dat.Amp         = Amp;
    dat.Orig        = Orig;
    dat.Gsolid0     = Gs0;
    dat.Gsolid1     = Gs1;
    dat.rev         = rev;
    dat.revrad1     = revrad1;
    dat.revrad2     = revrad2;



    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        
        dat.Bottom.Push(Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->ID);
        dat.Top   .Push(Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->ID);
    }
    //for (size_t i=0;i<ny;i++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    //}

    // Set inner drop
    int obsX = nx/2, obsY = 3.0*ny/4.0;
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

