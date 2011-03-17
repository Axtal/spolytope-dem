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
    Vec3_t          force;  // Force on planes
    double              S;  // Vertical separation of the planes
    std::ofstream  oss_ss;  // File to store the forces
};

void Setup (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.force = 0.5*(dat.p2->F-dat.p1->F);
    dat.S     = dat.p2->x(1)-dat.p1->x(1);
}

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    String fv;
    fv.Printf    ("%s_%08d_particles.res",dom.FileKey.CStr(), dom.idx_out);
    std::ofstream FV(fv.CStr());
    FV <<  Util::_10_6 << "PID" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" <<"\n";

    for (size_t i=0; i<dom.Particles.Size(); i++)
    {
        if (dom.Particles[i]->IsFree())
        {
            FV << Util::_8s << dom.Particles[i]->Index << Util::_8s << dom.Particles[i]->v(0)   << Util::_8s << dom.Particles[i]->v(1) << Util::_8s << dom.Particles[i]->v(2)
               << Util::_8s << dom.Particles[i]->x(0)  << Util::_8s << dom.Particles[i]->x(1)   << Util::_8s << dom.Particles[i]->x(2) << "\n";
        }
    }
    FV.close();
    size_t Nb  = dom.BInteractons.Size(); //Number of bonds
    size_t Nbb = 0;                   //Number of broken bonds
    if (Nb>0)
    {
        String fc;
        fc.Printf    ("%s_%08d_bonds.res",dom.FileKey.CStr(), dom.idx_out);
        std::ofstream FC(fc.CStr());
        FC <<  Util::_10_6 << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" <<  Util::_8s << "Area" << Util::_8s << "valid" << Util::_8s << "P1" <<  Util::_8s << "P2" <<"\n";

        for (size_t i=0; i<dom.BInteractons.Size(); i++)
        {
            FC << Util::_8s << dom.BInteractons[i]->Fnet(0)   << Util::_8s << dom.BInteractons[i]->Fnet(1)   << Util::_8s <<  dom.BInteractons[i]->Fnet(2) 
               << Util::_8s << dom.BInteractons[i]->Area      << Util::_8s << dom.BInteractons[i]->valid
               << Util::_8s << dom.BInteractons[i]->P1->Index << Util::_8s << dom.BInteractons[i]->P2->Index << "\n";
            if (!dom.BInteractons[i]->valid) Nbb++; //Add one broken bond if it is no longer valid
        }
        FC.close();

        String fl;
        fl.Printf    ("%s_%08d_clusters.res",dom.FileKey.CStr(), dom.idx_out);
        std::ofstream FL(fl.CStr());
        FL <<  Util::_10_6 << "Cluster" << Util::_8s << "Mass" << "\n";
        FL <<  Util::_8s   << 0         << Util::_8s << dom.Ms << "\n";
        for (size_t i=0;i<dom.Listofclusters.Size();i++)
        {
            double m = 0.0; //mass of the cluster
            for (size_t j=0;j<dom.Listofclusters[i].Size();j++)
            {
                m+=dom.Particles[dom.Listofclusters[i][j]]->Props.m;
            }
            FL << Util::_8s << i+1 << Util::_8s << m << "\n"; 
        }
        FL.close();
    }
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz" << Util::_8s << "ez \n";
    }
    else 
    {
        if (!dom.Finished) 
        {
            dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << fabs(dat.force(0)) << Util::_8s << fabs(dat.force(1)) << Util::_8s << fabs(dat.force(2)) << Util::_8s << dat.S << std::endl;
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


    d.CamPos = Vec3_t(0.1*Lx, 0.7*(Lx+Ly+Lz), 0.15*Lz); // position of camera
    if (ptype=="voronoi")      d.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.1);
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        d.GenFromMesh (mesh,/*R*/R,/*rho*/rho,true,false);
    }
    else throw new Fatal("Packing for particle type not implemented yet");
    d.Alpha = 0.5*R;
    d.Center();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    d.AddPlane(-2, Vec3_t(0.0,0.0,Xmin(2)-R), R, Lx, Ly, rho);
    d.AddPlane(-3, Vec3_t(0.0,0.0,Xmax(2)+R), R, Lx, Ly, rho);
    
    // properties of particles prior the brazilian test
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,Mu );
    B.Set(-2,"Bn Bt Bm Gn Gt eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,0.0);
    B.Set(-3,"Bn Bt Bm Gn Gt eps Kn Kt Mu",Bn,Bt,Bm,Gn,Gt,eps,Kn,Kt,0.0);
    d.SetProps(B);

    Vec3_t velocity(0.0,0.0,strf*(Xmax(2)-Xmin(2)+2*R)/Tf);

    Particle * p1 = d.GetParticle(-2);
    Particle * p2 = d.GetParticle(-3);
    p1->FixVeloc();
    p1->v =  velocity;
    p2->FixVeloc();
    p2->v = -velocity;
    dat.p1=p1;
    dat.p2=p2;

    d.WriteBPY(filekey.CStr());

    d.Solve(Tf, dt, dtOut, &Setup, &Report, filekey.CStr());


    return 0;
}
MECHSYS_CATCH

