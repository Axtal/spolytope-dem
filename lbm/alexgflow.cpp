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
    Array<Cell *> zmin;
    Array<Cell *> zmax;
    double          Kn;
    Vec3_t           g;
    Vec3_t        Xmin;
    Vec3_t        Xmax;
    double          Tf;
    double          dt;
    double         DPz;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    //double dRhobc;
    //(dom.Time > dat.Tf/4.0)&&(dom.Time < dat.Tf/2.0) ? dRhobc = 4.0*dat.DPz*dat.dt/(dat.Tf) : dRhobc = 0.0;

    
    for (size_t i=0;i<dat.zmin.Size();i++)
    {
        Cell * c = dat.zmin[i];
        if(c->IsSolid) continue;
        //c->RhoBC += dRhobc;
        c->F[5] = 1/3.0*(-2*c->F[0] - 2*c->F[1] - 4*c->F[12] - 4*c->F[13] - 2*c->F[2] - 2*c->F[3] - 2*c->F[4] - c->F[6] -  4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[7] = 1/24.0*(-2*c->F[0] + c->F[1] - 4*c->F[12] - 4*c->F[13] - 5*c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] +  20*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] - 4*c->F[13] + c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] + 20*c->F[9] + 2*c->RhoBC);
        c->F[11] = 1/24.0*(-2*c->F[0] + c->F[1] + 20*c->F[12] - 4*c->F[13] - 5*c->F[2] - 5*c->F[3] + c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14] = 1/24.0*(-2*c->F[0] - 5*c->F[1] - 4*c->F[12] + 20*c->F[13] + c->F[2] + c->F[3] - 5*c->F[4] - 4*c->F[6] - 4*c->F[8] - 4*c->F[9] + 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
        //std::cout << c->RhoBC << " " << dRhobc << " " << c->Rho << " " << dom.Time << std::endl;
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

    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Ff = dom.Particles[i]->Props.m*dat.g;
        double delta;
        delta =   dat.Xmin(0) - dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(0) += dat.Kn*delta;
        delta = - dat.Xmax(0) + dom.Particles[i]->x(0) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(0) -= dat.Kn*delta;
        delta =   dat.Xmin(1) - dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) + dom.Particles[i]->x(1) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) -= dat.Kn*delta;
        delta =   dat.Xmin(2) - dom.Particles[i]->x(2) + dom.Particles[i]->Props.R;
        if (delta > 0.0)  dom.Particles[i]->Ff(2) += dat.Kn*delta;
    }
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("report.res");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "C1" << Util::_8s << "C2" << Util::_8s << "Cn1" << Util::_8s << "Porosity" << Util::_8s << "Zmin" << Util::_8s << "Zmax" << std::endl;
    }
    if (!dom.Finished)
    {
        Vec3_t Xmin,Xmax;
        dom.BoundingBox(Xmin,Xmax);
        double volumecontainer = (Xmax(0)-Xmin(0))*(Xmax(1)-Xmin(1))*(Xmax(2)-Xmin(2));
        double por             = 1.0 - dom.Voltot/volumecontainer;
        double C1              = 0.0;
        double C2              = 0.0;
        double Cn              = 0.0;
        for (size_t i=0;i<dom.CInteractons.Size();i++)
        {
            DEM::CInteracton * CI = dom.CInteractons[i];
            if ((CI->P1->Tag==-1)||(CI->P2->Tag==-1))            
            {
                C1+=norm(CI->Fnet);
            }
            if ((CI->P1->Tag==-2)||(CI->P2->Tag==-2))            
            {
                C2+=norm(CI->Fnet);
            }
        }
        size_t n1 = 0;
        size_t n2 = 0;
        for(size_t i=0;i<dom.Particles.Size();i++)
        {
            Cn += dom.Particles[i]->Cn;
            if (dom.Particles[i]->Tag==-1) n1++;
            if (dom.Particles[i]->Tag==-2) n2++;
        }        
        C1/=n1;
        C2/=n2;
        Cn/=(n1+n2);
        dat.oss_ss << dom.Time << Util::_8s << C1 << Util::_8s << C2 << Util::_8s << Cn << Util::_8s << por << Util::_8s << Xmin(2) << Util::_8s << Xmax(2) << std::endl;
    }
    else
    {
        dat.oss_ss.close();
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

    bool   Render   = true;
    size_t seed     = 100;
    double fraction = 0.8;
    double Rmin     = 0.9;
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
    double R1       = 2.0;
    double R2       = 20.0;

    infile >> Render;    infile.ignore(200,'\n');
    infile >> seed;      infile.ignore(200,'\n');
    infile >> fraction;  infile.ignore(200,'\n');
    infile >> Rmin;      infile.ignore(200,'\n');
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
    infile >> R1;        infile.ignore(200,'\n');
    infile >> R2;        infile.ignore(200,'\n');



    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Step     = 2;
    dat.g        = 0.0,0.0,-g;
    dat.Xmin     = 0.0,0.0,dx*nx/5.0;
    dat.Xmax     = nx*dx,ny*dx,nz*dx;
    dat.Kn       = Kn;
    dat.Tf       = Tf;
    dat.dt       = dt;
    dat.DPz      = DPz;

    Dom.Alpha    = std::min(R1,R2);
    
    //Set solid boundaries
    //for (size_t i=0;i<nx;i++)
    //for (size_t j=0;j<nz;j++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
    //}
    //for (size_t i=0;i<ny;i++)
    //for (size_t j=0;j<nz;j++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j))->IsSolid = true;
    //}
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        dat.zmin.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0  )));
        dat.zmax.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,nz-1)));
    }
    
    //Initializing values
    double rho0 = 1.0;
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, OrthoSys::O);
    }

    //Setting boundary conditions
    for (size_t i=0;i<dat.zmin.Size();i++)
    {
        //if(dat.zmin[i]->Index(0)>nx/3.0&&dat.zmin[i]->Index(0)<2.0*nx/3.0&&dat.zmin[i]->Index(1)>ny/3.0&&dat.zmin[i]->Index(1)<2.0*ny/3.0)
        //{
            //dat.zmin[i]->RhoBC = rho0 + DPz;
            //dat.zmin[i]->Initialize(rho0 + DPz,OrthoSys::O);
        //}
        //else
        //{
            //dat.zmin[i]->IsSolid = true;
            //dat.zmin[i]->RhoBC = rho0;
            //dat.zmin[i]->Initialize(rho0,OrthoSys::O);
        //}
        //if(dat.zmax[i]->Index(0)>nx/3.0&&dat.zmax[i]->Index(0)<2.0*nx/3.0&&dat.zmax[i]->Index(1)>ny/3.0&&dat.zmax[i]->Index(1)<2.0*ny/3.0)
        //{
            //dat.zmax[i]->RhoBC = rho0;
            //dat.zmax[i]->Initialize(rho0,OrthoSys::O);
        //}
        //else
        //{
            //dat.zmax[i]->IsSolid = true;
            //dat.zmax[i]->RhoBC = rho0;
            //dat.zmax[i]->Initialize(rho0,OrthoSys::O);
        //}
        dat.zmin[i]->RhoBC = rho0 + DPz;
        dat.zmin[i]->Initialize(rho0 + DPz,OrthoSys::O);
        dat.zmax[i]->RhoBC = rho0;
        dat.zmax[i]->Initialize(rho0,OrthoSys::O);
    }
    {
        DEM::Domain DemDom; 
        DemDom.Load(filekey.CStr());
        
        Array<size_t> Big;
        Array<size_t> Small;
        for (size_t i=0;i<DemDom.Particles.Size();i++)
        {
            DEM::Particle * Pa = DemDom.Particles[i];
            if (Pa->Props.R>0.5*R2)  Big.Push(i);
            else                     Small.Push(i);
        }

        std::cout << DemDom.Particles.Size() << " " << Big.Size() << " " << Small.Size() << std::endl;
        size_t nmin = 0;
        size_t nmax = 0;
        size_t ratio = Small.Size()/Big.Size();
        size_t j;
        bool valid = true;
        for (size_t i=0;i<DemDom.Particles.Size();i++)
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
            DEM::Particle * Pa = DemDom.Particles[j];
            Dom.AddSphere(Pa->Tag,Pa->x,Pa->Props.R,2.5);
        }
    }

    //Dom.GenSpheresBox (1, dat.Xmin , dat.Xmax - Vec3_t(0.0,0.0,0.5*nz*dx), /*R*/R1, 2.5, seed, fraction, Rmin); ///< Create an array of spheres
    //Dom.GenSpheresBox (-2, dat.Xmax - Vec3_t(nx*dx,ny*dx,0.5*nz*dx), dat.Xmax + Vec3_t(0.0,0.0,0.3*nz*dx), /*R*/R2, 2.5, seed, fraction, Rmin); ///< Create an array of spheres
    
    
    for (size_t i=0;i<Dom.Particles.Size();i++)
    {
        Dom.Particles[i]->Props.Kn  =     Kn;
        Dom.Particles[i]->Props.Kt  =     Kn;
        Dom.Particles[i]->Props.Mu  =     Mu;
        Dom.Particles[i]->Props.Eta =    Eta;
        Dom.Particles[i]->Props.Beta=   Beta;
        Dom.Particles[i]->Props.Gn  =     Gn;
        Dom.Particles[i]->Props.Gt  =    0.0;
    }

    //Dom.WriteXDMF(filekey.CStr());

    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render,Nproc);
    return 0;
}
MECHSYS_CATCH


