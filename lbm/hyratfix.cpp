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
    double          Head;       ///< Current hydraulic head
    double            Tf;
    double            Kn;
    double           ome;
    double         dtOut;
    double          time;
    Vec3_t          Xmin;
    Vec3_t          Xmax;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup(Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat.Cells.Size();i++)
    {
        Cell * c = dom.Lat.Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
    if (dom.Time>dat.time)
    {
        dat.time += dat.dtOut;
    }
    //double rho = dat.Head*fabs(sin(dat.ome*dom.Time))+54.0;
    double rho = dat.Head*fabs(sin(dat.ome*dat.time))+54.0;
    for (size_t i=0;i<dat.Bottom.Size();i++)
    {
        dat.Bottom[i]->Initialize(rho,OrthoSys::O);
    }

}

void Report(Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double water = 0.0;
    double Sr    = 0.0;
    for (size_t i=0;i<dom.Lat.Cells.Size();i++)
    {
        Cell * c = dom.Lat.Cells[i];
        water   += c->Density();
        if (c->Density()>500.0) Sr+=1.0;
    }
    Sr/=(dom.Lat.Cells.Size()*(1-dom.Lat.SolidFraction()));
    double head = 0.0;
    for (size_t i=0;i<dom.Lat.Ndim(0);i++)
    {
        head += dom.Lat.GetCell(iVec3_t(i,1,0))->Density();
    }
    head/=dom.Lat.Ndim(0);
    double rho = dat.Head*fabs(sin(dat.ome*dom.Time))+54.0;
    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << head << Util::_8s << water << Util::_8s << Sr << std::endl;
}

int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());

    size_t nx   = 200;
    size_t ny   = 200;
    double Gs   =-200;
    double nu   = 0.01;
    double dx   = 1.0;
    double dt   = 1.0;
    double Tf   = 10000.0;
    double dtOut= 50.0;
    double rho  = 200.0;
    double ome  = 2.0;
    double Head = 500.0;

    {
        infile >> nx;    infile.ignore(200,'\n');
        infile >> ny;    infile.ignore(200,'\n');
        infile >> Gs;    infile.ignore(200,'\n');
        infile >> nu;    infile.ignore(200,'\n');
        infile >> dx;    infile.ignore(200,'\n');
        infile >> dt;    infile.ignore(200,'\n');
        infile >> Tf;    infile.ignore(200,'\n');
        infile >> dtOut; infile.ignore(200,'\n');
        infile >> rho;   infile.ignore(200,'\n'); 
        infile >> ome;   infile.ignore(200,'\n');
        infile >> Head;  infile.ignore(200,'\n');
    }

    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Lat.G    = -200.0;
    Dom.Lat.Gs   =  Gs;
    dat.g        = 0.0,-0.001,0.0;
    dat.Tf       = Tf;
    dat.Xmin     = 0.0,0.0,0.0;
    dat.Xmax     = nx*dx,ny*dx,0.0;
    dat.Kn       = 1.0e4*rho/500.0;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.Head     = Head;
    dat.dtOut    = dtOut;
    dat.time     = 0.0;

    //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    {
        dat.Bottom.Push(Dom.Lat.GetCell(iVec3_t(i,0   ,0)));
        Dom.Lat.GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat.GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat.GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(  54.0,v0);
        //if (j<ny/8.0) Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(1300.0,v0);
        //else          Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize( 54.0,v0);
    }

	// Set grains
	Table grains;
	grains.Read("circles.out");
	for (size_t i=0; i<grains["Xc"].Size(); ++i)
	{
        double DX = nx*dx/2.0;
        for (size_t j=0;j<2;j++)
        {
            for (size_t k=0;k<2;k++)
            {
		        double xc = grains["Xc"][i]*DX+j*DX;
		        double yc = grains["Yc"][i]*DX+(k+1)*DX;
		        double r  = grains["R" ][i]*DX*0.9;
                Dom.Lat.SolidDisk(Vec3_t(xc,yc,0.0),r);
            }
        }
		//double xc = grains["Xc"][i]*nx*dx;
		//double yc = grains["Yc"][i]*nx*dx+2.0*DX;
		//double r  = grains["R" ][i]*nx*dx*0.9;
        //Dom.Lat.SolidDisk(Vec3_t(xc,yc,0.0),r);
	}

    //Solving
    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(),std::ios::out);
    dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "PDen" << Util::_8s << "Head" << Util::_8s << "Water" << Util::_8s << "Sr" <<std::endl;
    Dom.Solve(Tf,dtOut,Setup,Report,"hyratfix");
    dat.oss_ss.close();
}
MECHSYS_CATCH

