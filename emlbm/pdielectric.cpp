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
#include <mechsys/dem/domain.h>
#include <mechsys/emlbm2/Domain.h>

using std::cout;
using std::endl;
struct UserData
{
    double ome;
    double Tf;
    double J0;
    int nx;
    int ny;
    int nz;
};

void Setup(EMLBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(dom.Nproc)
#endif
    for (int i=0;i<dat.nx;i++)
    for (int j=0;j<dat.ny;j++)
    for (int k=0;k<dat.nz;k++)
    {
        double dx=i;
        //dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e2*dat.J0*sin(dat.ome*dom.Time)*exp(-0.75*(dx*dx+dy*dy))*(tanh(double(k)-2.0*dat.nz/5.0)+tanh(3.0*dat.nz/5.0-double(k)));
        dom.Lat.GetCell(iVec3_t(i,j,k))->Jf = OrthoSys::e2*dat.J0*sin(dat.ome*dom.Time)*exp(-0.75*(dx*dx));
    }
}

void Report (EMLBM::Domain & dom, void *UD)
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

    String fileDEM;
    bool   Render   = true;
    size_t N        = 200;
    double ESol     = 5.0;
    double Tf       = 10000.0;
    double dtOut    = 50.0;
    double J0       = 1.0e-4;
    double ome      = 2.0;
    int    outlimit = 1;
    size_t buffer   = 1;
    {
        infile >> fileDEM;   infile.ignore(200,'\n');
        infile >> Render;    infile.ignore(200,'\n');
        infile >> N;         infile.ignore(200,'\n');
        infile >> ESol;      infile.ignore(200,'\n');
        infile >> Tf;        infile.ignore(200,'\n');
        infile >> dtOut;     infile.ignore(200,'\n');
        infile >> J0;        infile.ignore(200,'\n');
        infile >> ome;       infile.ignore(200,'\n');
        infile >> outlimit;  infile.ignore(200,'\n'); 
        infile >> buffer;    infile.ignore(200,'\n'); 
    }

    DEM::Domain DemDom;
    DemDom.Load(fileDEM.CStr());
    Array<int> idx(6);
    idx = -2,-3,-4,-5,-6,-7;
    DemDom.DelParticles(idx);
    Vec3_t Xmin,Xmax;
    DemDom.BoundingBox(Xmin,Xmax);
    //int    bound = -2;
    int    bound = outlimit;
    double dx = (Xmax(0)-Xmin(0))/(N-2*bound);
    size_t Ny = (Xmax(1)-Xmin(1))/dx + 2*bound;
    size_t Nz = (Xmax(2)-Xmin(2))/dx + 2*bound;
    DemDom.Center(0.5*(Xmax-Xmin)+Vec3_t(bound*dx,bound*dx,bound*dx));
    EMLBM::Domain Dom(iVec3_t(N,Ny,Nz), 1.0, 1.0);
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<Ny;j++)
        {
            for (int k=0;k<Nz;k++)
            {
                Vec3_t pos((i+0.5)*dx,(j+0.5)*dx,(k+0.5)*dx);
                for (size_t n=0;n<DemDom.Particles.Size();n++)
                {
                    DEM::Particle *P = DemDom.Particles[n];
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat.GetCell(iVec3_t(i,j,k))->Eps = ESol;
                    }
                }
            }
        }
    }
    Dom.Step = buffer;

    UserData dat;
    Dom.UserData = &dat;

    dat.Tf       = Tf;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.J0       = J0;
    dat.nx = N;
    dat.ny = Ny;
    dat.nz = Nz;

    //Dom.WriteXDMF("test");

    //Solving
    Dom.Solve(Tf,dtOut,&Setup,Report,filekey.CStr(),Render,Nproc);
}
MECHSYS_CATCH

