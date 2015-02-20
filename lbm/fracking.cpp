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
#include <mechsys/lbm/Domain.h>

using std::cout;
using std::endl;
struct UserData
{
    double            Orig;
    double             Amp;
    double              Tf;
    Array<Cell *>   Center;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double rho;
    rho = dat.Orig + dom.Time*(dat.Amp - dat.Orig)/dat.Tf;
    for(size_t i=0;i<dat.Center.Size();i++)
    {
        dat.Center[i]->Initialize(rho,OrthoSys::O);
    }
}

void Report (LBM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (!dom.Finished) 
    {
        String ff;
        ff.Printf    ("%s_bf_%04d",  dom.FileKey.CStr(), dom.idx_out);
        dom.WriteBF  (ff.CStr());
        ff.Printf    ("%s_frac_%04d",dom.FileKey.CStr(), dom.idx_out);
        dom.WriteFrac(ff.CStr());
    }
}

int main(int argc, char **argv) try
{
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    std::ifstream infile(filename.CStr());
    size_t Nproc = 1;
    if (argc>=3) Nproc = atoi(argv[2]);

    double Amax  = 10.0;
    size_t nx    = 500;
    size_t ny    = 500;
    size_t nz    = 500;
    double nu    = 0.2;
    double Orig  = 1.0;
    double Amp   = 1.0;
    double Tf    = 5000.0;
    double dt    = 1.0;
    double Alpha = 10.0;
    double sx    = 10.0;          
    double sy    = 20.0;          
    double r2    = 5.0;
    double Kn    = 3000.0;
    double Kt    = 3000.0;
    double Gn    = -0.2;
    double Mu    = 0.3;
    double Bn    = 10.0;
    double Bt    = 10.0;
    double Eps   = 0.01;

    infile >>  Amax ;       infile.ignore(200,'\n');
    infile >>  nx   ;       infile.ignore(200,'\n');
    infile >>  ny   ;       infile.ignore(200,'\n');
    infile >>  nz   ;       infile.ignore(200,'\n');
    infile >>  nu   ;       infile.ignore(200,'\n');
    infile >>  Orig ;       infile.ignore(200,'\n');
    infile >>  Amp  ;       infile.ignore(200,'\n');
    infile >>  Tf   ;       infile.ignore(200,'\n');
    infile >>  dt   ;       infile.ignore(200,'\n');
    infile >>  Alpha;       infile.ignore(200,'\n');
    infile >>  sx   ;       infile.ignore(200,'\n');
    infile >>  sy   ;       infile.ignore(200,'\n');
    infile >>  r2   ;       infile.ignore(200,'\n');
    infile >>  Kn   ;       infile.ignore(200,'\n');
    infile >>  Kt   ;       infile.ignore(200,'\n');
    infile >>  Gn   ;       infile.ignore(200,'\n');
    infile >>  Mu   ;       infile.ignore(200,'\n');
    infile >>  Bn   ;       infile.ignore(200,'\n');
    infile >>  Bt   ;       infile.ignore(200,'\n');
    infile >>  Eps  ;       infile.ignore(200,'\n');

    double r1 = 0.5*r2;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), 1.0, dt);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Dilate   = true;
    dat.Tf          = Tf;
    dat.Amp         = Amp;
    dat.Orig        = Orig;
    Dom.Alpha       = Alpha;

    size_t n_divisions = 10;

    Mesh::Unstructured mesh(/*NDim*/2);
    mesh.Set    (4+n_divisions,4+n_divisions, 1, 1);                            // 8 points, 8 segments, 1 region, 1 hole
    mesh.SetReg (0, -1, Amax, 0.2, 0.2);                 // id, tag, max{area}, x, y <<<<<<< regions
    mesh.SetHol (0, 0.5*nx, 0.5*ny);                     // id, x, y <<<<<<< holes
    mesh.SetPnt (0, -1, 0.0, 0.0);                       // id, vtag, x, y <<<<<< points
    mesh.SetPnt (1, -2,  nx, 0.0);                       // id, vtag, x, y
    mesh.SetPnt (2, -3,  nx,  ny);                       // id, vtag, x, y
    mesh.SetPnt (3, -4, 0.0,  ny);                       // id, vtag, x, y
    mesh.SetSeg (0, -10,  0, 1);                         // id, etag, L, R <<<<<<<<<<<< segments
    mesh.SetSeg (1, -20,  1, 2);                         // id, etag, L, R
    mesh.SetSeg (2, -30,  2, 3);                         // id, etag, L, R
    mesh.SetSeg (3, -40,  3, 0);                         // id, etag, L, R
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetPnt(i+4, 0, r2*cos(2*i*M_PI/n_divisions)+0.5*nx,r2*sin(2*i*M_PI/n_divisions)+0.5*ny);
    }
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetSeg(i+4, 0, i + 4, (i+1)%n_divisions + 4);
    }

    mesh.Generate ();
    Dom.GenFromMesh(mesh,/*spheroradius*/Alpha,/*density*/3.0,/*iscohesive*/true,/*montecarlo mass properties*/false,/*thickness*/double(nz));
    Dom.GenBoundingBox (/*InitialTag*/-2, 0.2*Alpha, /*Cf*/1.2, /*Rho*/3.0);
    Dom.Center(Vec3_t(0.5*(nx-1),0.5*(ny-1),0.5*(nz-1)));
    Array<int> DeleteTags(2);
    DeleteTags[0]  = -6;
    DeleteTags[1]  = -7;
    Dom.DelParticles(DeleteTags);
	
    Dict B;
    B.Set(-1,"Bn Bt Gn Eps Kn Kt Mu"   ,Bn,Bt,Gn,Eps,Kn,Kt,Mu );
    B.Set(-2,"Bn Bt Gn Eps Kn Kt Mu"   ,Bn,Bt,Gn,Eps,Kn,Kt,0.0);
    B.Set(-3,"Bn Bt Gn Eps Kn Kt Mu"   ,Bn,Bt,Gn,Eps,Kn,Kt,0.0);
    B.Set(-4,"Bn Bt Gn Eps Kn Kt Mu"   ,Bn,Bt,Gn,Eps,Kn,Kt,0.0);
    B.Set(-5,"Bn Bt Gn Eps Kn Kt Mu"   ,Bn,Bt,Gn,Eps,Kn,Kt,0.0);
    Dom.SetProps(B);

    Dom.GetParticle(-2)->Ff = -double(ny*nz)*sx*OrthoSys::e0;
    Dom.GetParticle(-3)->Ff =  double(ny*nz)*sx*OrthoSys::e0;
    Dom.GetParticle(-4)->Ff = -double(nx*nz)*sy*OrthoSys::e1;
    Dom.GetParticle(-5)->Ff =  double(nx*nz)*sy*OrthoSys::e1;
    Dom.GetParticle(-2)->FixVeloc();
    Dom.GetParticle(-3)->FixVeloc();
    Dom.GetParticle(-4)->FixVeloc();
    Dom.GetParticle(-5)->FixVeloc();
    Dom.GetParticle(-2)->vxf = false;
    Dom.GetParticle(-3)->vxf = false;
    Dom.GetParticle(-4)->vyf = false;
    Dom.GetParticle(-5)->vyf = false;

    
    for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
	for (size_t k=0; k<nz; ++k)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(Orig,OrthoSys::O);
		if (pow((int)(i)-0.5*nx,2.0) + pow((int)(j)-0.5*ny,2.0) <= pow(r1,2.0)) // circle equation
		{
            dat.Center.Push(Dom.Lat[0].GetCell(iVec3_t(i,j,k)));
		}
    }

    Dom.Solve(Tf,0.01*Tf,&Setup,&Report,"fracking",true,Nproc);


    return 0;
}
MECHSYS_CATCH

