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
// Sinking disks

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>

struct UserData
{
    Vec3_t             g;
    Array<Cell *> Bottom;
    Array<Cell *>    Top;
    double          Head;       ///< Current hydraulic head
    double          Orig;       ///< Original hydraulic head
    double            Tf;
    double            Kn;
    double           ome;
    double         dtOut;
    double          time;
    Vec3_t          Xmin;
    Vec3_t          Xmax;
    std::ofstream oss_ss;       ///< file for stress strain data
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
    if (dom.Time>dat.time)
    {
        dat.time += dat.dtOut;
    }
    double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
    for (size_t i=0;i<dat.Bottom.Size();i++)
    {
        dat.Bottom[i]->Initialize(rho,OrthoSys::O);
        //Cell * c = dat.Bottom[i];
		//if (c->IsSolid) continue;
        //c->RhoBC = rho;
		//double vy = -1.0 + (c->F[0]+c->F[1]+c->F[3] + 2.0*(c->F[4]+c->F[7]+c->F[8]))/c->RhoBC;
		//c->F[2] = c->F[4] - (2.0/3.0)*c->RhoBC*vy; 
		//c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vy - 0.5*(c->F[3]-c->F[1]);
		//c->F[5] = c->F[7] - (1.0/6.0)*c->RhoBC*vy + 0.5*(c->F[3]-c->F[1]);
        //c->Rho = c->VelDen(c->Vel);
    }
	// Cells with prescribed density
	for (size_t i=0; i<dat.Top.Size(); ++i)
	{
		Cell * c = dat.Top[i];
		if (c->IsSolid) continue;
		double vy = -1.0 + (c->F[0]+c->F[1]+c->F[3] + 2.0*(c->F[2]+c->F[5]+c->F[6]))/c->RhoBC;
		c->F[4] = c->F[2] - (2.0/3.0)*c->RhoBC*vy; 
		c->F[8] = c->F[6] - (1.0/6.0)*c->RhoBC*vy + 0.5*(c->F[3]-c->F[1]);
		c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vy - 0.5*(c->F[3]-c->F[1]);
        c->Rho = c->VelDen(c->Vel);
    }
}

void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double water = 0.0;
    double Sr    = 0.0;
    size_t nw    = 0;
    for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        if (c->Rho>1119.185) 
        {
            Sr+=1.0;
            water+=c->Rho;
            nw++;
        }
    }
    Sr/=(dom.Lat[0].Cells.Size()*(1-dom.Lat[0].SolidFraction()));
    if (nw>0) water/=nw;
    double head = 0.0;
    double Gasf = 0.0;
    double hmax = 0.0;
    for (size_t i=0;i<dom.Lat[0].Ndim(0);i++)
    {
        head += dom.Lat[0].GetCell(iVec3_t(i,1,0))->Rho;
        Gasf += dat.Top[i]->Rho*dat.Top[i]->Vel(1);
        double hmaxp = 0.0;
        for (size_t j=0;j<dom.Lat[0].Ndim(1);j++)
        {
            Cell * c = dom.Lat[0].GetCell(iVec3_t(i,j,0));
            if (c->Rho>1119.185&&j>hmaxp) hmaxp = j;
        }
        hmax += hmaxp;
    }
    head/=dom.Lat[0].Ndim(0);
    Gasf/=dom.Lat[0].Ndim(0);
    hmax/=dom.Lat[0].Ndim(0);
    double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << head << Util::_8s << water << Util::_8s << Sr << Util::_8s << hmax << Util::_8s << Gasf << std::endl;
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
    size_t nx       = 200;
    size_t ny       = 200;
    double Gs       = -200;
    double nu       = 0.05;
    double dx       = 1.0;
    double dt       = 1.0;
    double Tf       = 10000.0;
    double dtOut    = 50.0;
    double HeadStep = 1000.0;
    double rho      = 200.0;
    double ome      = 2.0;
    double Head     = 500.0;
    double Orig     = 54.0;
    double por      = 0.6;
    double seed     = 1000;

    {
        infile >> Render;    infile.ignore(200,'\n');
        infile >> nx;        infile.ignore(200,'\n');
        infile >> ny;        infile.ignore(200,'\n');
        infile >> Gs;        infile.ignore(200,'\n');
        infile >> nu;        infile.ignore(200,'\n');
        infile >> dx;        infile.ignore(200,'\n');
        infile >> dt;        infile.ignore(200,'\n');
        infile >> Tf;        infile.ignore(200,'\n');
        infile >> dtOut;     infile.ignore(200,'\n');
        infile >> HeadStep;  infile.ignore(200,'\n');
        infile >> rho;       infile.ignore(200,'\n'); 
        infile >> ome;       infile.ignore(200,'\n');
        infile >> Head;      infile.ignore(200,'\n');
        infile >> Orig;      infile.ignore(200,'\n');
        infile >> por;       infile.ignore(200,'\n');
        infile >> seed;      infile.ignore(200,'\n');
    }
    Array<double> nua(2);
    nua[0] = nu;
    nua[1] = nu;

    LBM::Domain Dom(D2Q9, nua, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Lat[0].G    = -200.0;
    Dom.Lat[0].Gs   =  Gs;
    dat.g        = 0.0,-0.001,0.0;
    dat.Tf       = Tf;
    dat.Xmin     = 0.0,0.0,0.0;
    dat.Xmax     = nx*dx,ny*dx,0.0;
    dat.Kn       = 1.0e4*rho/500.0;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.Head     = Head;
    dat.Orig     = Orig;
    dat.dtOut    = HeadStep;
    dat.time     = 0.0;

    Dom.Lat[1].G = 0.0;
    Dom.Lat[1].Gs= 0.0;
    Dom.Gmix     = 0.001;


    //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    {
        dat.Bottom.Push(Dom.Lat[0].GetCell(iVec3_t(i,0   ,0)));
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->Gs      = 0.0;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->Gs      = 0.0;
        dat.Top.Push   (Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0)));
        dat.Top[i]->RhoBC = rho;
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
    double        rc = 0.0;
    while (1-Dom.Lat[0].SolidFraction()/n>por)
    {
        ntries++;
        if (ntries>1.0e4) throw new Fatal("Too many tries to achieved requested porosity, please increase it");
        double Rmax = 0.02;
        double Rmin = 0.5*Rmax;
        double r  = ((Rmin*Rmax/(Rmax - double(rand())/RAND_MAX*(Rmax - Rmin))))*nx*dx;
        double DY = 0.0;
        double yc = DY + (ny*dx/6.0 - DY)*double(rand())/RAND_MAX+r+dx;
        double xc = nx*dx*double(rand())/RAND_MAX;
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

    ntries = 0;
    n      = 1.0;
    Radii.Resize(0);
    Xs.Resize(0);

    while (1-Dom.Lat[0].SolidFraction()/n>por)
    {
        ntries++;
        if (ntries>1.0e4) throw new Fatal("Too many tries to achieved requested porosity, please increase it");
        double Rmax = 0.04;
        double Rmin = 0.3*Rmax;
        double r  = ((Rmin*Rmax/(Rmax - double(rand())/RAND_MAX*(Rmax - Rmin))))*nx*dx;
        double DY = 1.0/6.0*ny*dx;
        double yc = DY + (ny*dx - DY)*double(rand())/RAND_MAX;
        double xc = nx*dx*double(rand())/RAND_MAX;
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
    std::cout << 1-Dom.Lat[0].SolidFraction() << std::endl;

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(  0.1,v0);
        Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(  rho,v0);
    }

    //Solving
    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(),std::ios::out);
    dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "PDen" << Util::_8s << "Head" << Util::_8s << "Water" << Util::_8s << "Sr" << Util::_8s << "Hmax" << Util::_8s << "Gf" << std::endl;
    Dom.Solve(Tf,dtOut,Setup,Report,"hyratfix",Render,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH

