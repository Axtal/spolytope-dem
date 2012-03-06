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
// Vorticity test

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/dem/domain.h>

struct UserData
{
    Array<Cell *> xmin0;
    Array<Cell *> xmax0;
    Array<Cell *> ymin0;
    Array<Cell *> ymax0;
    Array<Cell *> zmin0;
    Array<Cell *> zmax0;
    Array<Cell *> xmin1;
    Array<Cell *> xmax1;
    Array<Cell *> ymin1;
    Array<Cell *> ymax1;
    Array<Cell *> zmin1;
    Array<Cell *> zmax1;
    double         Head;       ///< Current hydraulic head
    double         Orig;       ///< Original hydraulic head
    double           Tf;
    double           Kn;
    double          ome;
    double        dtOut;
    double         time;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.Time>dat.time)
    {
        dat.time += dat.dtOut;
    }
    double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
    for (size_t i=0;i<dat.xmin0.Size();i++)
    {
        dat.xmin0[i]->Initialize(rho,OrthoSys::O);
    }
    for (size_t i=0;i<dat.xmax1.Size();i++)
    {
        Cell * c = dat.xmax1[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double water = 0.0;
    double oil   = 0.0;
    double Sr    = 0.0;
    size_t nw    = 0;
    size_t no    = 0;
    for (size_t i=0;i<dom.Lat[1].Cells.Size();i++)
    {
        Cell * c = dom.Lat[1].Cells[i];
        if (c->Rho>1119.185) 
        {
            Sr+=1.0;
            water+=c->Rho;
            nw++;
        }
        c = dom.Lat[0].Cells[i];
        if (c->Rho>800.0)
        {
            oil  +=c->Rho;
            no++;
        }
    }
    Sr/=(dom.Lat[0].Cells.Size()*(1-dom.Lat[0].SolidFraction()));
    if (nw>0) water/=nw;
    if (no>0) oil  /=no;
    double head = 0.0;
    size_t nfb  = 0;
    for (size_t i=0;i<dom.Lat[0].Ndim(1);i++)
    for (size_t j=0;j<dom.Lat[0].Ndim(2);j++)
    {
        Cell * c = dom.Lat[0].GetCell(iVec3_t(1,i,j));
        if (c->IsSolid) continue;
        head += c->Rho;        
        nfb++;
    }
    head/=nfb;
    double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << head << Util::_8s << water << Util::_8s << Sr << std::endl;
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
    size_t N        = 200;
    double Gs       = -200;
    double nu       = 0.05;
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
        infile >> N;         infile.ignore(200,'\n');
        infile >> Gs;        infile.ignore(200,'\n');
        infile >> nu;        infile.ignore(200,'\n');
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

    DEM::Domain DemDom;
    //DemDom.AddVoroPack(-1,R,10.0,10.0,10.0,nx,ny,nz,1.0,true,false,seed,fraction,Vec3_t(q,q,q));
    DemDom.GenSpheres  (-1, 10.0, 10, 1.0, "Normal", seed, 1.0);
    DemDom.Initialize (1.0);
    




    Vec3_t Xmin,Xmax;
    DemDom.BoundingBox(Xmin,Xmax);
    int    bound = -1;
    double dx = (Xmax(0)-Xmin(0))/(N-2*bound);
    double dy = (Xmax(1)-Xmin(1))/(N-2*bound);
    double dz = (Xmax(2)-Xmin(2))/(N-2*bound);
    DemDom.Center(0.5*(Xmax-Xmin)+Vec3_t(bound*dx,bound*dy,bound*dz));

    LBM::Domain Dom(D3Q15, nua, iVec3_t(N,N,N), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    dat.Tf       = Tf;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.Head     = Head;
    dat.Orig     = Orig;
    dat.dtOut    = HeadStep;
    dat.time     = 0.0;

    Dom.Lat[0].G = -150.0;
    Dom.Lat[0].Gs= -100.0;
    Dom.Lat[1].G = -200.0;
    Dom.Lat[1].Gs= -200.0;
    Dom.Gmix     =  0.0001;

	// set inner obstacle
    for (int i=0;i<N;i++)
    {
        for (int j=0;j<N;j++)
        {
            dat.xmin0.Push(Dom.Lat[0].GetCell(iVec3_t(0  ,i,j)));
            dat.xmax0.Push(Dom.Lat[0].GetCell(iVec3_t(N-1,i,j)));
            dat.ymin0.Push(Dom.Lat[0].GetCell(iVec3_t(i,0  ,j)));
            dat.ymax0.Push(Dom.Lat[0].GetCell(iVec3_t(i,N-1,j)));
            dat.zmin0.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0  )));
            dat.zmax0.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,N-1)));
            dat.xmin1.Push(Dom.Lat[1].GetCell(iVec3_t(0  ,i,j)));
            dat.xmax1.Push(Dom.Lat[1].GetCell(iVec3_t(N-1,i,j)));
            dat.ymin1.Push(Dom.Lat[1].GetCell(iVec3_t(i,0  ,j)));
            dat.ymax1.Push(Dom.Lat[1].GetCell(iVec3_t(i,N-1,j)));
            dat.zmin1.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,0  )));
            dat.zmax1.Push(Dom.Lat[1].GetCell(iVec3_t(i,j,N-1)));

            dat.xmax0.Last()->IsSolid = true;
            dat.xmin1.Last()->IsSolid = true;

            dat.ymin0.Last()->IsSolid = true;
            dat.ymax0.Last()->IsSolid = true;
            dat.zmin0.Last()->IsSolid = true;
            dat.zmax0.Last()->IsSolid = true;
            dat.ymin1.Last()->IsSolid = true;
            dat.ymax1.Last()->IsSolid = true;
            dat.zmin1.Last()->IsSolid = true;
            dat.zmax1.Last()->IsSolid = true;

            for (int k=0;k<N;k++)
            {
                Vec3_t pos((i+0.5)*dx,(j+0.5)*dy,(k+0.5)*dz);
                for (size_t n=0;n<DemDom.Particles.Size();n++)
                {
                    if (DemDom.Particles[n]->IsInsideAlt(pos))
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                }
            }
        }
    }

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(0.1, OrthoSys::O);
        Dom.Lat[1].Cells[i]->Initialize(1300.0, OrthoSys::O);
    }


    //Setting boundary conditions
    for (size_t i=0;i<dat.xmin0.Size();i++)
    {
        dat.xmax1[i]->RhoBC = 1300.0;
    }

    //Solving
    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(),std::ios::out);
    dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "PDen" << Util::_8s << "Head" << Util::_8s << "Water" << Util::_8s << "Sr" << std::endl;
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH


