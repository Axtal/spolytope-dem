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
    Particle         * p1;  // Upper plane
    Particle         * p2;  // Lower plane
    Array<Vec3_t *>    vm;  // pointer to the selected vertices
    Vec3_t             L0;
    Vec3_t          force;  // Force on planes
    std::ofstream  oss_ss;  // File to store the forces
};

void Setup (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.force = 0.5*(dat.p2->F-dat.p1->F);
}

void ReportSimple(DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "ez" << Util::_8s << "nu" << Util::_8s << "nu_ave" <<std::endl;
        double tol = dat.L0(0)/40.0;
        for (size_t i=0;i<dom.Particles.Size();i++)
        {
            for (size_t j=0;j<dom.Particles[i]->Verts.Size();j++)
            {
                double r = sqrt(pow((*dom.Particles[i]->Verts[j])(0),2) + pow((*dom.Particles[i]->Verts[j])(2),2));
                if (fabs(r - 0.5*dat.L0(0))<tol)
                {
                    dat.vm.Push(dom.Particles[i]->Verts[j]);
                }
            }
        }
    }
    else 
    {
        if (!dom.Finished) 
        {
            Vec3_t Xmin, Xmax;
            dom.BoundingBox(Xmin, Xmax);
            double S      = (Xmax(1) - Xmin(1) - dat.L0(1))/dat.L0(1);
            std::cout << Xmin << Xmax << dat.L0 << std::endl;
            double nu     = -0.5*((Xmax(0) - Xmin(0) - dat.L0(0))/dat.L0(0) + (Xmax(2) - Xmin(2) - dat.L0(2))/dat.L0(2))/S;
            double nu_ave = 0.0;
            for (size_t i=0;i<dat.vm.Size();i++)
            {
                double r = sqrt(pow((*dat.vm[i])(0),2) + pow((*dat.vm[i])(2),2));
                nu_ave += 2.0*r;
            }
            nu_ave = -(nu_ave/dat.vm.Size() - dat.L0(0))/(dat.L0(0)*S);
            dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << fabs(dat.force(0)) << Util::_8s << fabs(dat.force(1)) << Util::_8s << fabs(dat.force(2)) << Util::_8s << S << Util::_8s << nu << Util::_8s << nu_ave <<  std::endl;
        }
        else
        {
            dat.oss_ss.close();
        }
    }
}

int main(int argc, char **argv) try
{
    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    ifstream infile(filename.CStr());
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    UserData dat; 
    Domain d(&dat);
    String ptype;       // Particle type 
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Bn          = 1.0e6;
    double Bt          = 3.3e5;
    double Bm          = 3.3e5;
    double eps         = 0.01;
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt          = 0.0001;
    double dtOut       = 0.1;
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double Tf          = 10.0;
    double strf        = 0.01;

    infile >> ptype;              infile.ignore(200,'\n');
    infile >> Kn;                 infile.ignore(200,'\n');
    infile >> Kt;                 infile.ignore(200,'\n');
    infile >> Gn;                 infile.ignore(200,'\n');
    infile >> Gt;                 infile.ignore(200,'\n');
    infile >> Mu;                 infile.ignore(200,'\n');
    infile >> Bn;                 infile.ignore(200,'\n');
    infile >> Bt;                 infile.ignore(200,'\n');
    infile >> Bm;                 infile.ignore(200,'\n');
    infile >> eps;                infile.ignore(200,'\n');
    infile >> R;                  infile.ignore(200,'\n');
    infile >> seed;               infile.ignore(200,'\n');
    infile >> dt;                 infile.ignore(200,'\n');
    infile >> dtOut;              infile.ignore(200,'\n');
    infile >> Lx;                 infile.ignore(200,'\n');
    infile >> Ly;                 infile.ignore(200,'\n');
    infile >> Lz;                 infile.ignore(200,'\n');
    infile >> nx;                 infile.ignore(200,'\n');
    infile >> ny;                 infile.ignore(200,'\n');
    infile >> nz;                 infile.ignore(200,'\n');
    infile >> rho;                infile.ignore(200,'\n');
    infile >> Tf;                 infile.ignore(200,'\n');
    infile >> strf;               infile.ignore(200,'\n');


    double cam_x=1.1*Ly, cam_y=0.0, cam_z=0.0;
    d.CamPos = cam_x, cam_y, cam_z;
    if (ptype=="voronoi")      d.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.1);
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        //mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        size_t n_divisions = 30;
        mesh.Set(2*n_divisions,n_divisions+2, 1, 0);
        mesh.SetReg (0,  -1, 0.25*M_PI*Lx*Lx*Ly/(0.3*nx*ny*nz),  0.0, 0.0, 0.0);  // id, tag, max{volume}, x, y, z <<<<<<< regions
        //mesh.SetReg (0,  -1, 2.0,  0.0, 0.0, 0.0);  // id, tag, max{volume}, x, y, z <<<<<<< regions
        Array<int> F1;
        Array<int> F2;
        for(size_t i=0; i<n_divisions; i++)
        {
            mesh.SetPnt(i              , 0, 0.5*Lx*cos(2*i*M_PI/n_divisions), -0.5*Ly, 0.5*Lx*sin(2*i*M_PI/n_divisions));
            mesh.SetPnt(i+n_divisions  , 0, 0.5*Lx*cos(2*i*M_PI/n_divisions),  0.5*Ly, 0.5*Lx*sin(2*i*M_PI/n_divisions));
            F1.Push(i);
            F2.Push(2*n_divisions-1-i);
        }
        mesh.SetFac(0,0,F1);
        mesh.SetFac(1,0,F2);
        for(size_t i=0; i<n_divisions; i++)
        {
            //mesh.SetSeg(i               , 0, i            , (i+1)%n_divisions);
            //mesh.SetSeg(i+   n_divisions, 0, i+n_divisions, (i+1)%n_divisions+n_divisions);
            //mesh.SetSeg(i+ 2*n_divisions, 0, i            ,  i+n_divisions);
            Array<int> F(4);
            F = i+n_divisions, (i+1)%n_divisions + n_divisions, (i+1)%n_divisions, i;
            mesh.SetFac(i+2, 0, F);
        }
        mesh.Generate();
        d.GenFromMesh (mesh,/*R*/R,/*rho*/rho,true,false);
    }
    else throw new Fatal("Packing for particle type not implemented yet");
    d.Alpha = 0.5*R;
    d.Center();
    d.GenBoundingPlane(-2,R,1.0,true);
    
    // properties of particles prior the brazilian test
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt ,     eps,Kn,Kt,Mu );
    B.Set(-2,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,0.0,-0.1*eps,Kn,Kt,0.0);
    B.Set(-3,"Bn Bt Bm Gn Gt Eps Kn Kt Mu",Bn,Bt,Bm,Gn,0.0,-0.1*eps,Kn,Kt,0.0);
    d.SetProps(B);


    Particle * p1 = d.GetParticle(-2);
    Particle * p2 = d.GetParticle(-3);
    p1->FixVeloc();
    p2->FixVeloc();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    Vec3_t velocity(0.0,strf*(Xmax(1)-Xmin(1))/Tf,0.0);
    dat.L0 = Xmax-Xmin;
    p1->v =  0.5*velocity;
    p2->v = -0.5*velocity;
    dat.p1=p1;
    dat.p2=p2;

    d.WriteBPY(filekey.CStr());

    d.Solve(Tf, dt, dtOut, &Setup, &ReportSimple, filekey.CStr());


    return 0;
}
MECHSYS_CATCH

