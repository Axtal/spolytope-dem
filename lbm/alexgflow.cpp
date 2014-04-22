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

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/lbm/Domain.h>

struct UserData
{
    Array<Cell *> zmin;
    Array<Cell *> zmax;
    double          Kn;
    Vec3_t           g;
    double          Tf;
    double          dt;
    double         DPz;
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    for (size_t i=0;i<dat.zmin.Size();i++)
    {
        Cell * c = dat.zmin[i];
        if(c->IsSolid) continue;
        c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
        c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.zmax.Size();i++)
    {
        Cell * c = dat.zmax[i];
        if(c->IsSolid) continue;
        c->F[6]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]- 2*c->F[2]- 2*c->F[3]- 2*c->F[4]- c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[14]+ c->F[2]- 5*c->F[3]+ c->F[4]- 4*c->F[5]+ 20*c->F[7]+ 2*c->RhoBC);
        c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[14]- 5*c->F[2]+ c->F[3]- 5*c->F[4]- 4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[14]+ c->F[2]+ c->F[3]- 5*c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[14]- 5*c->F[2]- 5*c->F[3]+ c->F[4]-4*c->F[5]- 4*c->F[7]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }

}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc==3) Nproc = atoi(argv[2]);

    String fileDEM;
    bool   Render   = true;
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
    double Mu       = 0.4;
    double Eta      = 1.0;
    double Beta     = 0.12;
    double Gn       = 1.0;
    double DPz      = 0.0;

    infile >> fileDEM;   infile.ignore(200,'\n');
    infile >> Render;    infile.ignore(200,'\n');
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
    infile >> Mu;        infile.ignore(200,'\n');
    infile >> Eta;       infile.ignore(200,'\n');
    infile >> Beta;      infile.ignore(200,'\n');
    infile >> Gn;        infile.ignore(200,'\n');
    infile >> DPz;       infile.ignore(200,'\n');


    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Step     = 2;
    dat.g        = 0.0,0.0,-g;
    dat.Kn       = Kn;
    dat.Tf       = Tf;
    dat.dt       = dt;
    dat.DPz      = DPz;

    
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        dat.zmin.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0  )));
        dat.zmax.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1)));
    }
    
    //Initializing values
    double rho0 = 1000.0;
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        rho0 = 1000.0 + (DPz*(nz-Dom.Lat[0].Cells[i]->Index(2)))/nz;
        Dom.Lat[0].Cells[i]->Initialize(rho0, OrthoSys::O);
    }

    //Setting boundary conditions
    for (size_t i=0;i<dat.zmin.Size();i++)
    {
        dat.zmin[i]->RhoBC = rho0 + DPz;
        dat.zmin[i]->Initialize(rho0 + DPz,OrthoSys::O);
        dat.zmax[i]->RhoBC = rho0;
        dat.zmax[i]->Initialize(rho0,OrthoSys::O);
    }

    //Loading DEM sample
    Dom.Load(fileDEM.CStr());

    //Setting fixed walls that do not interact with the fluid
    Dom.GetParticle(-3)->FixVeloc();
    Dom.GetParticle(-4)->FixVeloc();
    Dom.GetParticle(-5)->FixVeloc();
    Dom.GetParticle(-6)->FixVeloc();
    Dom.GetParticle(-7)->FixVeloc();
    Dom.GetParticle(-3)->Bdry = true;
    Dom.GetParticle(-4)->Bdry = true;
    Dom.GetParticle(-5)->Bdry = true;
    Dom.GetParticle(-6)->Bdry = true;
    Dom.GetParticle(-7)->Bdry = true;

    Vec3_t Xmin,Xmax;
    Dom.BoundingBox(Xmin,Xmax);
    Vec3_t Cen = 0.5*(Xmax-Xmin)+10.0*dx*OrthoSys::e2;
    Dom.Center(Cen);


    double maxR = 0.0; 
    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        Dom.Particles[i]->Props.Kn  =     Kn;
        Dom.Particles[i]->Props.Kt  =     Kn;
        Dom.Particles[i]->Props.Mu  =     Mu;
        Dom.Particles[i]->Props.Eta =    Eta;
        Dom.Particles[i]->Props.Beta=   Beta;
        Dom.Particles[i]->Props.Gn  =     Gn;
        Dom.Particles[i]->Props.Gt  =    0.0;
        maxR = std::max(Dom.Particles[i]->Props.R,maxR);
    }

    Array<size_t> Big;
    Array<size_t> Small;
    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        DEM::Particle * Pa = Dom.Particles[i];
        if (Pa->Props.R>0.5*maxR)  Big.Push(i);
        else                       Small.Push(i);
    }

    std::cout << Dom.Particles.Size() << " " << Big.Size() << " " << Small.Size() << std::endl;
    size_t nmin = 0;
    size_t nmax = 0;
    size_t ratio = Small.Size()/Big.Size();
    size_t j;
    bool valid = true;
    Array<DEM::Particle * >  TempPar(Dom.Particles.Size());
    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        if (nmin%ratio==0&&valid&&nmax<Big.Size()) 
        {
            j = Big[nmax];
            nmax++;
            valid = false;
        }
        else
        {
            j = Small[nmin];
            nmin++;
            valid = true;
        }
        DEM::Particle * Pa = Dom.Particles[j];
        TempPar[i] = Pa; 
    }

    Dom.Particles.Resize(0);
    Dom.Particles = TempPar;

    Dom.WriteXDMF(filekey.CStr());

    //Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    return 0;
}
MECHSYS_CATCH


