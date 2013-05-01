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

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;
using std::ofstream;
using DEM::Domain;

struct UserData
{
};

void Setup (DEM::Domain & dom, void *UD)
{
}

void Report (DEM::Domain & dom, void *UD)
{
}

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc>2) Nproc = atoi(argv[2]);

    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    UserData dat; 
    Domain d(&dat);
    size_t   Render;
    double Lx = 10.0;
    double Ly = 5.0;
    double theta = 45.0;
    double Amax  = 0.5;
    double thickness = 1.0;
    bool   Cohesion = true;
    double Kn          = 1.0e6;
    double Kt          = 3.3e5;
    double Bn          = 1.0e6;
    double Bt          = 3.3e5;
    double Bm          = 3.3e5;
    double Gn          = 16.0;
    double Gt          = 8.0;
    double Mu          = 0.3;
    double eps         = 0.01;
    double g           = 9.8;
    double Tf          = 10.0;
    double dt          = 1.0e-5;

    infile >> Render;          infile.ignore(200,'\n');
    infile >> Lx;              infile.ignore(200,'\n');
    infile >> Ly;              infile.ignore(200,'\n');
    infile >> theta;           infile.ignore(200,'\n');
    infile >> Amax;            infile.ignore(200,'\n');
    infile >> thickness;       infile.ignore(200,'\n');
    infile >> Cohesion;        infile.ignore(200,'\n');
    infile >> Kn ;             infile.ignore(200,'\n');
    infile >> Kt ;             infile.ignore(200,'\n');
    infile >> Bn ;             infile.ignore(200,'\n');
    infile >> Bt ;             infile.ignore(200,'\n');
    infile >> Bm ;             infile.ignore(200,'\n');
    infile >> Gn ;             infile.ignore(200,'\n');
    infile >> Gt ;             infile.ignore(200,'\n');
    infile >> Mu ;             infile.ignore(200,'\n');
    infile >> eps;             infile.ignore(200,'\n');
    infile >> g  ;             infile.ignore(200,'\n');
    infile >> Tf ;             infile.ignore(200,'\n');
    infile >> dt ;             infile.ignore(200,'\n');



    d.CamPos = Vec3_t(0.0, 0.0, 1.5*Lx); // position of camera

    Mesh::Unstructured mesh(2);
    double sbase = 0.5*Ly/tan(M_PI*theta/180.0);
    mesh.Set    (6, 6, 1, 0);  // division, edges, regions and holes

    mesh.SetReg (0,  -1,  Amax,0.5*Lx,0.25*Ly);  // id, tag, max{volume}, x, y, z <<<<<<< regions

    mesh.SetPnt(0,0,0.0,   0.0);
    mesh.SetPnt(1,0,Lx ,   0.0);
    mesh.SetPnt(2,0,Lx ,0.33*Ly);
    mesh.SetPnt(3,0,0.5*Lx,0.33*Ly);
    mesh.SetPnt(4,0,0.5*Lx-sbase,Ly);
    mesh.SetPnt(5,0,0.0,Ly);

    for (size_t i=0;i<6;i++)
    {
        mesh.SetSeg(i, 0, i, (i+1)%6);
    }

    mesh.Generate();
    d.GenFromMesh(mesh,0.1*sqrt(Amax/10),3.0,Cohesion,false,thickness);
    d.Alpha = 0.5*sqrt(Amax/10);
    d.Center();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    d.AddPlane(-2, Vec3_t(0.0,Xmin(1)-0.5*sqrt(Amax/10),0.0), 0.5*sqrt(Amax/10), 1.1*Lx, 1.2*thickness, 1.0, M_PI/2.0, &OrthoSys::e0);
    d.AddPlane(-3, Vec3_t(0.0,0.0,Xmin(2)-0.5*sqrt(Amax/10)), 0.5*sqrt(Amax/10), 1.1*Lx, 1.1*Ly, 1.0);
    d.AddPlane(-4, Vec3_t(0.0,0.0,Xmax(2)+0.5*sqrt(Amax/10)), 0.5*sqrt(Amax/10), 1.1*Lx, 1.1*Ly, 1.0);
    d.AddPlane(-5, Vec3_t(Xmin(0)-0.5*sqrt(Amax/10),0.0,0.0), 0.5*sqrt(Amax/10), 1.2*thickness, 1.1*Ly, 1.0, M_PI/2.0, &OrthoSys::e1);
    d.AddPlane(-6, Vec3_t(Xmax(0)+0.5*sqrt(Amax/10),0.0,0.0), 0.5*sqrt(Amax/10), 1.2*thickness, 1.1*Ly, 1.0, M_PI/2.0, &OrthoSys::e1);


    d.GetParticle(-2,true)->FixVeloc();
    d.GetParticle(-3,true)->FixVeloc();
    d.GetParticle(-4,true)->FixVeloc();
    d.GetParticle(-5,true)->FixVeloc();
    d.GetParticle(-6,true)->FixVeloc();
    
    // properties of particles prior the brazilian test
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,Mu);
    B.Set(-2,"Gn Gt Eps Kn Kt Mu",Gn,0.0,eps,Kn,Kt,0.0);
    B.Set(-3,"Gn Gt Eps Kn Kt Mu",Gn,0.0,eps,Kn,Kt,0.0);
    B.Set(-4,"Gn Gt Eps Kn Kt Mu",Gn,0.0,eps,Kn,Kt,0.0);
    B.Set(-5,"Gn Gt Eps Kn Kt Mu",Gn,0.0,eps,Kn,Kt,0.0);
    B.Set(-6,"Gn Gt Eps Kn Kt Mu",Gn,0.0,eps,Kn,Kt,0.0);
    d.SetProps(B);

    size_t n_div = 20;
    double step = Lx/n_div;

    for (size_t i=0;i<d.Particles.Size();i++)
    {
        if (!d.Particles[i]->IsFree()) continue;
        int tag = d.Particles[i]->x(0)/step + 0.5*n_div + 1;
        tag = tag%2;
        d.Particles[i]->Tag = tag;
        d.Particles[i]->Ff  = d.Particles[i]->Props.m*Vec3_t(0.0,-g,0.0);
    }
    
    d.WriteXDMF("slope");



    d.Solve(Tf, dt, 0.01*Tf, &Setup, &Report, filekey.CStr(),Render,Nproc);


    return 0;
}
MECHSYS_CATCH
