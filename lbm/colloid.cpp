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
// Colloid transport

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>

struct UserData
{
    Array<Cell *> Left; 
    Array<Cell *> Right;
    size_t         alt;
    size_t         Npr;
    size_t         Npl;
    size_t          nr;
    double        vmax; 
    double          Kn; 
    double        xlim;
    double          rc;
    double          vp;
    Vec3_t        Xmin; 
    Vec3_t        Xmax; 
    std::ofstream oss_ss;       ///< file for stress strain data
};

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

	// Cells with prescribed velocity
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
		if (c->IsSolid) continue;
		double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-c->VelBC(0));
		c->F[1] = c->F[3] + (2.0/3.0)*rho*c->VelBC(0);
		c->F[5] = c->F[7] + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC(1) - 0.5*(c->F[2]-c->F[4]);
		c->F[8] = c->F[6] + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC(1) + 0.5*(c->F[2]-c->F[4]);
	}

	// Cells with prescribed density
	for (size_t i=0; i<dat.Right.Size(); ++i)
	{
		Cell * c = dat.Right[i];
		if (c->IsSolid) continue;
		double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
	}

    dat.Npr = 0;
    dat.Npl = 0;
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        if (dom.Particles[i]->X(0) > dat.Xmax(0) && dom.Particles[i]->IsFree()) dat.Npr++;
        if (dom.Particles[i]->X(0) > dat.xlim    && dom.Particles[i]->IsFree()) dat.Npl++;
        dom.Particles[i]->Ff = 0.0,0.0,0.0;
        double delta;
        delta =   dat.Xmin(0) - dom.Particles[i]->X(0) + dom.Particles[i]->R;
        if (delta > 0.0)  dom.Particles[i]->Ff(0) += dat.Kn*delta;
        //delta = - dat.Xmax(0) + dom.Particles[i]->X(0) + dom.Particles[i]->R;
        //if (delta > 0.0)  dom.Particles[i]->Ff(0) -= dat.Kn*delta;
        delta =   dat.Xmin(1) - dom.Particles[i]->X(1) + dom.Particles[i]->R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) += dat.Kn*delta;
        delta = - dat.Xmax(1) + dom.Particles[i]->X(1) + dom.Particles[i]->R;
        if (delta > 0.0)  dom.Particles[i]->Ff(1) -= dat.Kn*delta;
    }

    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    if (Xmin(0) > dat.xlim + 2.0*dat.rc)
    {
        for (double y=0.05*dat.Xmax(1) + dat.alt*2.0*dat.rc;y<0.95*dat.Xmax(1);y+=4.0*dat.rc)
        {
            dom.AddDisk(0,Vec3_t(dat.xlim, y,0.0),Vec3_t(dat.vp,0.0,0.0),OrthoSys::O,3.0,dat.rc,1.0);
            dom.Particles[dom.Particles.Size()-1]->Kn = dat.Kn;
        }
        dom.ResetContacts();
        dom.ResetDisplacements();
        dat.alt ++;
        dat.alt = (dat.alt%2);
    }
}

void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("colloid_n_particles.res");
        dat.oss_ss.open(fs.CStr(),std::ios::out);
        dat.oss_ss << Util::_10_6  <<  "Time" << Util::_8s << "Npl" << Util::_8s << "Npr" << Util::_8s << "PGrad" << Util::_8s << "MassFlux" << Util::_8s << "Porosity" << std::endl;
    }
    double MassFlux = 0.0;
    double Density  = 0.0;
    for (size_t i = 0;i < dom.Lat[0].Ndim(1) ; i++)
    {
        Vec3_t V;
        double rho = dom.Lat[0].GetCell(iVec3_t(0,i,0))->VelDen(V);
        Density  += rho;
        MassFlux += rho*norm(V);
    }
    Density /=dom.Lat[0].Ndim(1);
    MassFlux/=dom.Lat[0].Ndim(1);

    double As = 0.0;
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        double r = dom.Particles[i]->R;
        if (dom.Particles[i]->X(0) - r > dat.nr*dat.rc) As += M_PI*r*r;
    }
    double porosity = 1.0 - As/((dat.Xmax(0) - dat.nr*dat.rc)*dat.Xmax(1));

    dat.oss_ss << dom.Time << Util::_8s << dat.Npl << Util::_8s << dat.Npr << Util::_8s << fabs(Density - 1.0) << Util::_8s << MassFlux << Util::_8s << porosity << std::endl;
}

int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());

    bool   Render = true;
    size_t seed   = 1000;
    size_t nx     = 400;
    size_t ny     = 400;
    double nu     = 0.1;
    double dx     = 1.0;
    double dt     = 1.0;    
    double vb     = 0.1;
    double por    = 0.5;
    double rc     = 5.0;
    double Kn     = 1.0e2;
    double Tf     = 1.0e6;
    double dtOut  = 1.0e4;
    {
        infile >> Render;    infile.ignore(200,'\n');
        infile >> seed  ;    infile.ignore(200,'\n');
        infile >> nx    ;    infile.ignore(200,'\n');
        infile >> ny    ;    infile.ignore(200,'\n');  
        infile >> nu    ;    infile.ignore(200,'\n');  
        infile >> dx    ;    infile.ignore(200,'\n');  
        infile >> dt    ;    infile.ignore(200,'\n');  
        infile >> vb    ;    infile.ignore(200,'\n');  
        infile >> por   ;    infile.ignore(200,'\n');  
        infile >> rc    ;    infile.ignore(200,'\n');  
        infile >> Kn    ;    infile.ignore(200,'\n');  
        infile >> Tf    ;    infile.ignore(200,'\n');  
        infile >> dtOut ;    infile.ignore(200,'\n');  
    }               


    double xlim   = 2.0*rc;
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Alpha    = 2.0*rc;
    dat.xlim     = xlim;
    dat.Xmin     = 0.0,0.0,0.0;
    dat.Xmax     = nx*dx,ny*dx,0.0;
    dat.rc       = rc;
    dat.vp       = 0.0*vb;
    dat.alt      = 1;
    dat.Kn       = Kn;
    dat.Npr      = 0;
    dat.Npl      = 0;
    dat.nr       = 10;

    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0)));
        double vx = vb; // horizontal velocity
        double vy = 0.0;                          // vertical velocity
		Vec3_t v(vx, vy, 0.0);                    // velocity vector
        dat.Left [i]->VelBC = v;
        dat.Left [i]->RhoBC = 1.0;
        dat.Right[i]->VelBC = v;
        dat.Right[i]->RhoBC = 1.0;

    }

    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }

	// Set grains
	//Table grains;
	//grains.Read("circles.out");
	//for (size_t i=0; i<grains["Xc"].Size(); ++i)
	//{
		//double xc = grains["Xc"][i]*nx+0.1*nx;
		//double yc = grains["Yc"][i]*ny;
		//double r  = grains["R" ][i]*nx*0.8;
        //Dom.AddDisk(0,Vec3_t(xc,yc,0.0),OrthoSys::O,OrthoSys::O,3.0,r,1.0);
        //Dom.Particles[Dom.Particles.Size()-1]->FixVelocity();
        //Dom.Particles[Dom.Particles.Size()-1]->Kn = Kn;
        //Dom.Particles[Dom.Particles.Size()-1]->ImprintDisk(Dom.Lat);
	//}
    srand(seed);
    size_t ntries = 0;
    double n      = (nx*dx - dat.nr*rc)/(nx*dx);
    while (1-Dom.Lat[0].SolidFraction()/n>por)
    {
        ntries++;
        if (ntries>1.0e4) throw new Fatal("Too many tries to achieved requested porosity, please increase it");
        double Rmin = 0.02;
        double Rmax = 0.10;
        double r  = ((Rmin*Rmax/(Rmax - double(rand())/RAND_MAX*(Rmax - Rmin))))*nx*dx;
        double DX = dat.nr*rc + r;
        double xc = DX + (nx*dx - DX)*double(rand())/RAND_MAX;
        double yc = 0.0*ny*dx + (ny*dx - 0.0*ny*dx)*double(rand())/RAND_MAX;
        Vec3_t X(xc,yc,0.0);
        bool invalid = false;
        for (size_t i=0;i<Dom.Particles.Size();i++)
        {
            //if (norm(Dom.Particles[i]->X - X) < std::min(r,Dom.Particles[i]->R))
            if (norm(Dom.Particles[i]->X - X) < r + Dom.Particles[i]->R + rc)
            {
                invalid = true;
                break;
            }
        }
        if (invalid) continue;
        Dom.AddDisk(0,Vec3_t(xc,yc,0.0),OrthoSys::O,OrthoSys::O,3.0,r,1.0);
        Dom.Particles[Dom.Particles.Size()-1]->FixVelocity();
        Dom.Particles[Dom.Particles.Size()-1]->Kn = Kn;
        Dom.ImprintLattice();
    }

    std::cout << 1-Dom.Lat[0].SolidFraction()/n << std::endl;

    //Writing correlation information
    String fs;
    fs.Printf("constriction.res");
    std::ofstream os;
    os.open(fs.CStr(),std::ios::out);
    os << Util::_10_6  << "CR" << Util::_8s << "X" << Util::_8s << "Y" <<std::endl;
    for (size_t i=0;i<Dom.Particles.Size()-1;i++)
    {
        for (size_t j=i+1;j<Dom.Particles.Size();j++)
        {
            double dist = norm(Dom.Particles[i]->X - Dom.Particles[j]->X);
            double r1   = Dom.Particles[i]->R;
            double r2   = Dom.Particles[j]->R;
            double cons = dist - r1 - r2;
            if (cons > 4*rc) continue;
            Vec3_t ri   = (r1 + 0.5*cons)*(Dom.Particles[j]->X - Dom.Particles[i]->X)/dist + Dom.Particles[i]->X;
            os << cons << Util::_8s << ri(0) << Util::_8s << ri(1) << std::endl;
        }
    }
    os.close();

    for (double y=0.05*dat.Xmax(1);y<0.95*dat.Xmax(1);y+=4.0*rc)
    {
        Dom.AddDisk(0,Vec3_t(xlim, y,0.0),Vec3_t(0.0,0.0,0.0),OrthoSys::O,3.0,rc,1.0);
        Dom.Particles[Dom.Particles.Size()-1]->Kn = Kn;
    }

    

    double rho0 = 1.0;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }

    //Solving
    Dom.Time = 0.0;

    Dom.Solve(Tf,dtOut,Setup,Report,"colloid",Render);
    dat.oss_ss.close();
 
}
MECHSYS_CATCH

