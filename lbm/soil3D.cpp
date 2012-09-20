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

//STD
#include<iostream>
#include <list>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/dem/domain.h>

struct UserData
{
    Array<Cell *>  xmin0;
    Array<Cell *>  xmax0;
    Array<Cell *>  xmin1;
    Array<Cell *>  xmax1;
    double          Head;       ///< Current hydraulic head
    double          Orig;       ///< Original hydraulic head
    double            Tf;
    double            Kn;
    double           ome;
    double         dtOut;
    double          time;
    double           rho;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.Orig>0.0)
    {
        if (dom.Time>dat.time)
        {
            dat.time += dat.dtOut;
        }
        //double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
        double a   = M_PI/dat.ome;
        double rho = dat.Head*((1.0/a)*(dat.time-a*(floor(dat.time/a)+0.5))*pow(-1.0,floor(dat.time/a))+0.5)+dat.Orig;
        double rho0min = 0.999*rho;
        double rho1min = 0.001*rho;
        double rho0max = 0.001*(2.0*dat.rho - rho);
        double rho1max = 0.999*(2.0*dat.rho - rho);
        for (size_t i=0;i<dat.xmin0.Size();i++)
        {
            Cell * c = dat.xmin0[i];
            c->RhoBC = rho0min;
            c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);

            c = dat.xmin1[i];
            c->RhoBC = rho1min;
            c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);

            c = dat.xmax0[i];
            c->RhoBC = rho0max;
            c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);
            
            c = dat.xmax1[i];
            c->RhoBC = rho1max;
            c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);
        }
    }
    else
    {
        if (dom.Time>=dat.time)
        {
            dat.time += M_PI/dat.ome;
            dat.Head *= -1.0;
        }
        for (size_t i=0;i<dat.xmin0.Size();i++)
        {
            Cell * c = dat.xmin0[i];
            c->RhoBC = 0.999*dat.Head + c->F[0]+2.0*c->F[2]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2.0*c->F[8]+2.0*c->F[10]+2.0*c->F[12]+2.0*c->F[14];
            c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);
            //std::cout << dom.Time << " " << i << " " << c->Rho << " " << c->RhoBC << " " << c->Vel << " " << dat.Head << " " << c->Rho*c->Vel(0) << std::endl;

            c = dat.xmin1[i];
            c->RhoBC = 0.001*dat.Head + c->F[0]+2.0*c->F[2]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2.0*c->F[8]+2.0*c->F[10]+2.0*c->F[12]+2.0*c->F[14];
            c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
            c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);

            c = dat.xmax0[i];
            c->RhoBC = 0.001*dat.Head + c->F[0]+2.0*c->F[1]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2.0*c->F[7]+2.0*c->F[9]+2.0*c->F[11]+2.0*c->F[13];
            c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);

            c = dat.xmax1[i];
            c->RhoBC = 0.999*dat.Head + c->F[0]+2.0*c->F[1]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2.0*c->F[7]+2.0*c->F[9]+2.0*c->F[11]+2.0*c->F[13];
            c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
            c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
            c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
            c->Rho = c->VelDen(c->Vel);
        }
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
        double wr = dom.Lat[1].Cells[i]->Rho;
        double ar = dom.Lat[0].Cells[i]->Rho;
        if (dom.Lat[1].Cells[i]->IsSolid) continue;
        if (wr>0.5*dat.rho) 
        {
            Sr+=1.0;
            water+=(wr + ar + dom.Gmix*wr*ar)/3.0;
            nw++;
        }
        if (ar>0.5*dat.rho)
        {
            oil  +=(wr + ar + dom.Gmix*wr*ar)/3.0;
            no++;
        }
    }
    Sr/=(dom.Lat[0].Cells.Size()*(1-dom.Lat[0].SolidFraction()));
    if (nw>0) water/=nw;
    if (no>0) oil  /=no;
    double rhow = 0.0;
    double rhoo = 0.0;
    size_t nfb  = 0;
    size_t nfo  = 0;
    for (size_t i=0;i<dom.Lat[0].Ndim(1);i++)
    for (size_t j=0;j<dom.Lat[0].Ndim(2);j++)
    {
        Cell * c = dom.Lat[0].GetCell(iVec3_t(1,i,j));
        if (c->IsSolid) continue;
        rhow += c->Rho;        
        nfb++;
        c = dom.Lat[1].GetCell(iVec3_t(dom.Lat[1].Ndim(0)-2,i,j));
        rhoo += c->Rho;        
        nfo++;
    }
    rhow/=nfb;
    rhoo/=nfo;
    double Pc; 
    double rho;
    
    if (dat.Orig>0.0)
    {
        //double rho = dat.Head*sin(dat.ome*dat.time)*sin(dat.ome*dat.time)+dat.Orig;
        double a   = M_PI/dat.ome;
        rho        = dat.Head*((1.0/a)*(dat.time-a*(floor(dat.time/a)+0.5))*pow(-1.0,floor(dat.time/a))+0.5)+dat.Orig;
        Pc         = (2.0*(rho - dat.rho) + dom.Gmix*(rho*rho*0.999*0.001 - (2.0*dat.rho - rho)*(2.0*dat.rho - rho)*0.999*0.001))/3.0;
    }
    else
    {
        double pxmin0 = 0.0;
        double pxmin1 = 0.0;
        double pxmax0 = 0.0;
        double pxmax1 = 0.0;
        for (size_t i=0;i<dat.xmin0.Size();i++)
        {
            Cell * c = dat.xmin0[i];
            pxmin0 += c->Rho;
            c = dat.xmin1[i];
            pxmin1 += c->Rho;
            c = dat.xmax0[i];
            pxmax0 += c->Rho;
            c = dat.xmax1[i];
            pxmax1 += c->Rho;
        }
        pxmin0/= dat.xmin0.Size();
        pxmin1/= dat.xmin0.Size();
        pxmax0/= dat.xmin0.Size();
        pxmax1/= dat.xmin0.Size();
        Pc         = (pxmin0 + pxmin1 - pxmax0 - pxmax1 + dom.Gmix*(pxmin0*pxmin1 - pxmax0*pxmax1))/3.0;
        rho        = pxmin0 + pxmin1;
    }

    dat.oss_ss << dom.Time << Util::_8s << rho << Util::_8s << rhoo << Util::_8s << rhow << Util::_8s << water << Util::_8s << oil << Util::_8s << Pc << Util::_8s << Sr << std::endl;
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
    String fileLBM;
    bool   Render   = true;
    size_t N        = 200;
    double Gs0      = -0.53;
    double Gs1      = -0.53;
    double Gmix     = 2.0;
    double nu       = 0.05;
    double dt       = 1.0;
    double Tf       = 10000.0;
    double dtOut    = 50.0;
    double HeadStep = 1000.0;
    double rho      = 200.0;
    double ome      = 2.0;
    double Head     = 500.0;
    double Orig     = 54.0;
    {
        infile >> fileDEM;   infile.ignore(200,'\n');
        infile >> fileLBM;   infile.ignore(200,'\n');
        infile >> Render;    infile.ignore(200,'\n');
        infile >> N;         infile.ignore(200,'\n');
        infile >> Gs0;       infile.ignore(200,'\n');
        infile >> Gs1;       infile.ignore(200,'\n');
        infile >> Gmix;      infile.ignore(200,'\n');
        infile >> nu;        infile.ignore(200,'\n');
        infile >> dt;        infile.ignore(200,'\n');
        infile >> Tf;        infile.ignore(200,'\n');
        infile >> dtOut;     infile.ignore(200,'\n');
        infile >> HeadStep;  infile.ignore(200,'\n');
        infile >> rho;       infile.ignore(200,'\n'); 
        infile >> ome;       infile.ignore(200,'\n');
        infile >> Head;      infile.ignore(200,'\n');
        infile >> Orig;      infile.ignore(200,'\n');
    }
    Array<double> nua(2);
    nua[0] = nu;
    nua[1] = nu;


    DEM::Domain DemDom;
    DemDom.Load(fileDEM.CStr());
    Array<int> idx(6);
    idx = -2,-3,-4,-5,-6,-7;
    DemDom.DelParticles(idx);
    Vec3_t Xmin,Xmax;
    DemDom.BoundingBox(Xmin,Xmax);
    int    bound = 4;
    double dx = (Xmax(0)-Xmin(0))/(N-2*bound);
    size_t Ny = (Xmax(1)-Xmin(1))/dx-bound;
    size_t Nz = (Xmax(2)-Xmin(2))/dx-bound;
    DemDom.Center(0.5*(Xmax-Xmin)+Vec3_t(bound*dx,0.0,0.0));
    LBM::Domain Dom(D3Q15, nua, iVec3_t(N,Ny,Nz), 1.0, 1.0);
    //Dom.PrtVec = false;
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
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                    Vec3_t T(0.0,-(Ny*dx),0.0);
                    P->Translate(T);
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                    T = Vec3_t(0.0,Ny*dx,-(Nz*dx));
                    P->Translate(T);
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                    T = Vec3_t(0.0,Ny*dx,Nz*dx);
                    P->Translate(T);
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                    T = Vec3_t(0.0,-(Ny*dx),Nz*dx);
                    P->Translate(T);
                    if (P->IsInsideAlt(pos)) 
                    {
                        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                        Dom.Lat[1].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    }
                    T = Vec3_t(0.0,0.0,-(Nz*dx));
                    P->Translate(T);
                }
            }
        }
    }

    UserData dat;
    Dom.UserData = &dat;

    dat.Tf       = Tf;
    dat.ome      = 2*M_PI*ome/Tf;
    dat.Orig     = Orig;
    dat.dtOut    = HeadStep;
    dat.time     = 0.0;
    dat.rho      = rho;
    if (dat.Orig<0.0)
    {
        //dat.Head = dat.Orig*dat.rho*(1.0-Dom.Lat[0].SolidFraction())*N*2.0*ome/Tf;
        dat.Head = dat.Orig;
    }
    else dat.Head     = Head;

    Dom.Lat[0].G = 0.0;
    Dom.Lat[0].Gs= Gs0;
    Dom.Lat[1].G = 0.0;
    Dom.Lat[1].Gs= Gs1;
    Dom.Gmix     = Gmix;

	// set inner obstacles
    
    for (int i=0;i<Ny;i++)
    for (int j=0;j<Nz;j++)
    {
        dat.xmin0.Push(Dom.Lat[0].GetCell(iVec3_t(0  ,i,j)));
        dat.xmax0.Push(Dom.Lat[0].GetCell(iVec3_t(N-1,i,j)));
        dat.xmin1.Push(Dom.Lat[1].GetCell(iVec3_t(0  ,i,j)));
        dat.xmax1.Push(Dom.Lat[1].GetCell(iVec3_t(N-1,i,j)));
    }

    
    //Initializing values
    if (Util::FileExists(fileLBM))
    {
        hid_t file_id;
        file_id = H5Fopen(fileLBM.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
        float * Density0   = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Vvec0      = new float[3*Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Density1   = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        float * Vvec1      = new float[3*Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
        H5LTread_dataset_float(file_id,"Density_0" ,Density0);
        H5LTread_dataset_float(file_id,"Velocity_0",Vvec0   );
        H5LTread_dataset_float(file_id,"Density_1" ,Density1);
        H5LTread_dataset_float(file_id,"Velocity_1",Vvec1   );
        for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
        {
            Cell * c = Dom.Lat[0].Cells[i];
            Vec3_t V;
            V(0) =  Vvec0[3*i  ];
            V(1) =  Vvec0[3*i+1];
            V(2) =  Vvec0[3*i+2];
            c->Initialize(Density0[i],V);
            c = Dom.Lat[1].Cells[i];
            V(0) =  Vvec1[3*i  ];
            V(1) =  Vvec1[3*i+1];
            V(2) =  Vvec1[3*i+2];
            c->Initialize(Density1[i],V);
        }
        delete [] Density0;
        delete [] Vvec0;
        delete [] Density1;
        delete [] Vvec1;
    }
    else
    {
        for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
        {
            Cell * c0 = Dom.Lat[0].Cells[i];
            Cell * c1 = Dom.Lat[1].Cells[i];
            if (c0->Index(0)>bound) 
            {
                c1->Initialize(0.999*rho, OrthoSys::O);
                c0->Initialize(0.001*rho, OrthoSys::O);
            }
            else
            {
                c0->Initialize(0.999*rho, OrthoSys::O);
                c1->Initialize(0.001*rho, OrthoSys::O);
            }
        }
    }

    //Solving
    String fs;
    fs.Printf("water_retention.res");
    dat.oss_ss.open(fs.CStr(),std::ios::out);
    dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "PDen" << Util::_8s << "Rhow" << Util::_8s << "Rhoo" << Util::_8s << "Water" << Util::_8s << "Oil" << Util::_8s << "Pc" << Util::_8s << "Sr" << std::endl;
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    dat.oss_ss.close();
}
MECHSYS_CATCH


