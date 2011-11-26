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
    Array<Cell *> xmin;
    Array<Cell *> xmax;
    Array<Cell *> ymin;
    Array<Cell *> ymax;
    Array<Cell *> zmin;
    Array<Cell *> zmax;
    double        vmax;
    double        rho;
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        Cell * c = dat.xmin[i];
        if(c->IsSolid) continue;
        c->F[1] = 1.0/3.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[7] = 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
        c->F[9] = 1.0/24.0*(-2*c->F[0]+20*c->F[10]-4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[11]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]+20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
        c->F[13]= 1.0/24.0*(-2*c->F[0]-4*c->F[10]-4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[1]  = (1.0/3.0)*(-c->F[2] + 2.0*(c->RhoBC - c->F[0] - c->F[3] - c->F[4] - c->F[5] - c->F[6]) - 4.0*(c->F[8] + c->F[10] + c->F[12] + c->F[14]));
        //c->F[7 ] = (1.0/24.0)*(c->F[3] + c->F[5] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[10] + c->F[12] + c->F[14]) - 5.0*(c->F[4] + c->F[6]) + 20.0*c->F[8 ]);
        //c->F[9 ] = (1.0/24.0)*(c->F[3] + c->F[6] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[12] + c->F[14]) - 5.0*(c->F[4] + c->F[5]) + 20.0*c->F[10]);
        //c->F[11] = (1.0/24.0)*(c->F[4] + c->F[5] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[10] + c->F[14]) - 5.0*(c->F[3] + c->F[6]) + 20.0*c->F[12]);
        //c->F[13] = (1.0/24.0)*(c->F[4] + c->F[6] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[10] + c->F[12]) - 5.0*(c->F[3] + c->F[5]) + 20.0*c->F[14]);
        //std::cout << c->Density() << std::endl;
        c->Rho = c->VelDen(c->Vel);
    }

    for (size_t i=0;i<dat.xmax.Size();i++)
    {
        Cell * c = dat.xmax[i];
        if(c->IsSolid) continue;
        c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
        c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
        c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        //c->F[2]  = (1.0/3.0)*(-c->F[1] + 2.0*(c->RhoBC - c->F[0] - c->F[3] - c->F[4] - c->F[5] - c->F[6]) - 4.0*(c->F[7] + c->F[9 ] + c->F[11] + c->F[13]));
        //c->F[8 ] = (1.0/24.0)*(c->F[3] + c->F[5] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[10] + c->F[12] + c->F[14]) - 5.0*(c->F[4] + c->F[6]) + 20.0*c->F[8 ]);
        //c->F[10] = (1.0/24.0)*(c->F[3] + c->F[6] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[12] + c->F[14]) - 5.0*(c->F[4] + c->F[5]) + 20.0*c->F[10]);
        //c->F[12] = (1.0/24.0)*(c->F[4] + c->F[5] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[10] + c->F[14]) - 5.0*(c->F[3] + c->F[6]) + 20.0*c->F[12]);
        //c->F[14] = (1.0/24.0)*(c->F[4] + c->F[6] + 2.0*(c->RhoBC - c->F[0]) - 4.0*(c->F[2] + c->F[8 ] + c->F[10] + c->F[12]) - 5.0*(c->F[3] + c->F[5]) + 20.0*c->F[14]);
        //std::cout << c->Density() << std::endl;
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.ymin.Size();i++)
    {
        Cell * c = dat.ymin[i];
        if(c->IsSolid) continue;
        c->F[3]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 2*c->F[2]- c->F[4]- 2*c->F[5]- 2*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[7]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]+ c->F[5]-  5*c->F[6]+ 20*c->F[8]+ 2*c->RhoBC);
        c->F[9]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[10]- 4*c->F[11]- 4*c->F[13]- 5*c->F[2]- 4*c->F[4]- 5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[12]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]+ 20*c->F[11]- 4*c->F[13]+ c->F[2]- 4*c->F[4]-5*c->F[5]+ c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->F[14]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[10]- 4*c->F[11]+ 20*c->F[13]+ c->F[2]- 4*c->F[4]+ c->F[5]-5*c->F[6]- 4*c->F[8]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
    for (size_t i=0;i<dat.ymax.Size();i++)
    {
        Cell * c = dat.ymax[i];
        if(c->IsSolid) continue;
        c->F[4]= 1/3.0*(-2*c->F[0]- 2*c->F[1]- 4*c->F[12]- 4*c->F[14]- 2*c->F[2]- c->F[3]- 2*c->F[5]- 2*c->F[6]- 4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[8]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]+  20*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[10]= 1/24.0*(-2*c->F[0]- 5*c->F[1]- 4*c->F[12]- 4*c->F[14]+ c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]- 4*c->F[7]+ 20*c->F[9]+ 2*c->RhoBC);
        c->F[11]= 1/24.0*(-2*c->F[0]+ c->F[1]+ 20*c->F[12]- 4*c->F[14]- 5*c->F[2]- 4*c->F[3]+ c->F[5]- 5*c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->F[13]= 1/24.0*(-2*c->F[0]+ c->F[1]- 4*c->F[12]+ 20*c->F[14]- 5*c->F[2]- 4*c->F[3]- 5*c->F[5]+ c->F[6]-4*c->F[7]- 4*c->F[9]+ 2*c->RhoBC);
        c->Rho = c->VelDen(c->Vel);
    }
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

    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

	// Cells with prescribed velocity
	//for (size_t i=0; i<dat.Left.Size(); ++i)
	//{
		//Cell * c = dat.Left[i];
		//if (c->IsSolid) continue;
		//double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-c->VelBC(0));
		//c->F[1] = c->F[3] + (2.0/3.0)*rho*c->VelBC(0);
		//c->F[5] = c->F[7] + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC[1] - 0.5*(c->F[2]-c->F[4]);
		//c->F[8] = c->F[6] + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC[1] + 0.5*(c->F[2]-c->F[4]);
	//}

	// Cells with prescribed density
	//for (size_t i=0; i<dat.Right.Size(); ++i)
	//{
		//Cell * c = dat.Right[i];
		//if (c->IsSolid) continue;
		//double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		//c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		//c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		//c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
	//}
}

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("flux.res");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "M" << std::endl;
    }
    Vec3_t Flux = OrthoSys::O;
    double M    = 0.0;
    size_t nc   = 0;
    for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
    {
        Cell * c = dom.Lat[0].Cells[i];
        if (c->IsSolid) continue;
        Vec3_t DF;
        double rho = c->VelDen(DF);
        Flux += rho*DF;
        M += rho;
        nc++;
    }
    Flux/=M;
    dat.oss_ss << dom.Time << Util::_8s << Flux(0) << Util::_8s << Flux(1) << Util::_8s << Flux(2) << Util::_8s << M/nc << std::endl;
}


int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());

    bool   Render = true;
    size_t N      = 100;
    size_t nx;
    size_t ny;
    size_t nz;
    size_t seed;
    double fraction;
    double q;
    double R;
    double nu     = 1.0/6.0;
    double Tf     = 1.0e3;
    double dtOut  = 1.0e1;
    bool   oct    = true;
    double Giso   = 1.0;
    double Gdev   = 3.0;
    double th     = 0.0;
    double DPx;
    double DPy;
    double DPz;

    infile >> Render;   infile.ignore(200,'\n');
    infile >> N;        infile.ignore(200,'\n');
    infile >> nx;       infile.ignore(200,'\n');
    infile >> ny;       infile.ignore(200,'\n');
    infile >> nz;       infile.ignore(200,'\n');
    infile >> seed;     infile.ignore(200,'\n');
    infile >> fraction; infile.ignore(200,'\n');
    infile >> q;        infile.ignore(200,'\n');
    infile >> R;        infile.ignore(200,'\n');
    infile >> nu;       infile.ignore(200,'\n');
    infile >> Tf;       infile.ignore(200,'\n');
    infile >> dtOut;    infile.ignore(200,'\n');
    infile >> oct;      infile.ignore(200,'\n');
    if (oct)
    {
        infile >> Giso;     infile.ignore(200,'\n');
        infile >> Gdev;     infile.ignore(200,'\n');
        infile >> th;       infile.ignore(200,'\n');
        DPx    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0-2.0*M_PI/3.0);
        DPy    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0);
        DPz    = Giso/3.0 + 2.0/3.0*Gdev*sin(M_PI*th/180.0+2.0*M_PI/3.0);
    }
    else
    {
        infile >> DPx;      infile.ignore(200,'\n');
        infile >> DPy;      infile.ignore(200,'\n');
        infile >> DPz;      infile.ignore(200,'\n');
    }
    

    DEM::Domain DemDom;
    DemDom.AddVoroPack(-1,R,10.0,10.0,10.0,nx,ny,nz,1.0,true,false,seed,fraction,Vec3_t(q,q,q));
    std::ofstream areafile("area.out");
    for (size_t i=0;i<DemDom.BInteractons.Size();i++)
    {
        Vec3_t Area;
        DemDom.BInteractons[i]->P1->Faces[DemDom.BInteractons[i]->F1]->Normal(Area);
        Area *= DemDom.BInteractons[i]->Area;
        areafile << Util::_8s << Area(0) << Util::_8s << Area(1) << Util::_8s << Area(2) << std::endl;
    }
    areafile.close();



    Vec3_t Inet(OrthoSys::O);
    for (size_t i=0;i<DemDom.Particles.Size();i++)
    {   
        list<double> II;
        II.push_back(DemDom.Particles[i]->I(0));
        II.push_back(DemDom.Particles[i]->I(1));
        II.push_back(DemDom.Particles[i]->I(2));
        II.sort();
        list<double>::iterator it;
        size_t j = 0;
        for (it = II.begin();it != II.end();it++)
        {
            Inet(j) += *it;
            j++;
        }
    }
    Inet/=DemDom.Particles.Size();
    std::ofstream parfile("param.inp");
    parfile << Util::_8s << DPx     << Util::_8s << DPy     << Util::_8s << DPz     << std::endl;
    parfile << Util::_8s << Inet(0) << Util::_8s << Inet(1) << Util::_8s << Inet(2) << std::endl;
    parfile << Util::_8s << nx      << Util::_8s << ny      << Util::_8s << nz      << std::endl;
    parfile.close();
    //Array<int> Tags(6);
    //Tags[0] = -2;
    //Tags[1] = -3;
    //Tags[2] = -4;
    //Tags[3] = -5;
    //Tags[4] = -6;
    //Tags[5] = -7;
    //DemDom.DelParticles(Tags);
    //double dx = Cf*lmin/N;
    for (size_t n=0;n<DemDom.Particles.Size();n++)
    {
        DemDom.Particles[n]->Props.R*=0.0;
    }
    Vec3_t Xmin,Xmax;
    DemDom.BoundingBox(Xmin,Xmax);
    int    bound = -1;
    double dx = (Xmax(0)-Xmin(0))/(N-2*bound);
    double dy = (Xmax(1)-Xmin(1))/(N-2*bound);
    double dz = (Xmax(2)-Xmin(2))/(N-2*bound);
    DemDom.Center(0.5*(Xmax-Xmin)+Vec3_t(bound*dx,bound*dy,bound*dz));
    //DemDom.BoundingBox(Xmin,Xmax);
    //std::cout << Xmin << Xmax << Vec3_t(dx,dy,dz)<< std::endl;
    
    
    LBM::Domain Dom(D3Q15, nu, iVec3_t(N,N,N), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

	// set inner obstacle
    for (int i=0;i<N;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0  ,0  ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,N-1,0  ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,0  ,N-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,N-1,N-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,i,0  ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,i,0  ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,i,N-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,i,N-1))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,0  ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,0  ,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(0  ,N-1,i))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(N-1,N-1,i))->IsSolid = true;
        for (int j=0;j<N;j++)
        {
            dat.xmin.Push(Dom.Lat[0].GetCell(iVec3_t(0  ,i,j)));
            dat.xmax.Push(Dom.Lat[0].GetCell(iVec3_t(N-1,i,j)));
            dat.ymin.Push(Dom.Lat[0].GetCell(iVec3_t(i,0  ,j)));
            dat.ymax.Push(Dom.Lat[0].GetCell(iVec3_t(i,N-1,j)));
            dat.zmin.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,0  )));
            dat.zmax.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,N-1)));
            for (int k=0;k<N;k++)
            {
                Vec3_t pos((i+0.5)*dx,(j+0.5)*dy,(k+0.5)*dz);
                for (size_t n=0;n<DemDom.Particles.Size();n++)
                {
                    if (DemDom.Particles[n]->IsInsideAlt(pos)) Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
                    //if (i==0&&j==75&&k==5) std::cout << DemDom.Particles[n]->IsInside(pos);
                }
            }
        }
    }

    //Initializing values
    double rho0 = 1.0;
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, OrthoSys::O);
    }

    //Setting boundary conditions
    for (size_t i=0;i<dat.xmin.Size();i++)
    {
        dat.xmin[i]->RhoBC = rho0 + DPx;
        dat.ymin[i]->RhoBC = rho0 + DPy;
        dat.zmin[i]->RhoBC = rho0 + DPz;
        dat.xmax[i]->RhoBC = rho0;
        dat.ymax[i]->RhoBC = rho0;
        dat.zmax[i]->RhoBC = rho0;
        dat.xmin[i]->Initialize(rho0+DPx,OrthoSys::O);
        dat.ymin[i]->Initialize(rho0+DPy,OrthoSys::O);
        dat.zmin[i]->Initialize(rho0+DPz,OrthoSys::O);
        dat.xmax[i]->Initialize(rho0,OrthoSys::O);
        dat.ymax[i]->Initialize(rho0,OrthoSys::O);
        dat.zmax[i]->Initialize(rho0,OrthoSys::O);
    }

    //Solving
    Dom.Solve(Tf,dtOut,Setup,Report,filekey.CStr(),Render);
    dat.oss_ss.close();
}
MECHSYS_CATCH


