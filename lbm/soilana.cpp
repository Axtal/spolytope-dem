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
    double ds       = 1.0;
    double st       = 100.0;
    double dthmin   = 0.8;
    double dthmax   = 0.9;
    size_t Nmin     = 0;
    size_t Nmax     = 0;
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
        infile >> ds;        infile.ignore(200,'\n');
        infile >> st;        infile.ignore(200,'\n');
        infile >> dthmin;    infile.ignore(200,'\n');
        infile >> dthmax;    infile.ignore(200,'\n');
        infile >> Nmin;      infile.ignore(200,'\n');
        infile >> Nmax;      infile.ignore(200,'\n');
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


    Dom.Lat[0].G = 0.0;
    Dom.Lat[0].Gs= Gs0;
    Dom.Lat[1].G = 0.0;
    Dom.Lat[1].Gs= Gs1;
    Dom.Gmix     = Gmix;

    Array<size_t> NoSolid;

    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Cell * c = Dom.Lat[0].Cells[i];
        bool nei_solid = false;
        if (c->IsSolid) nei_solid = true;
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb = Dom.Lat[0].Cells[c->Neighs[j]];
            if (nb->IsSolid) nei_solid = true;
            if (nb->ID>c->ID) 
            {
                if (!c->IsSolid||!nb->IsSolid) Dom.CellPairs.Push(iVec3_t(i,nb->ID,j));
            }
        }
        if (!nei_solid
            &&c->Index(0)>bound&&c->Index(0)<Dom.Lat[0].Ndim(0)-bound
            &&c->Index(1)>bound&&c->Index(1)<Dom.Lat[0].Ndim(1)-bound
            &&c->Index(2)>bound&&c->Index(2)<Dom.Lat[0].Ndim(2)-bound)
        {
            NoSolid.Push(c->ID);
        }
    }

    LBM::MtData MTD[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = &Dom;
        MTD[i].Dmx      = 0.0;
        MTD[i].dt       = Dom.Lat[0].dt;
    }
    pthread_t thrs[Nproc];   

    
    std::ofstream oss_ss;       ///< file for stress strain data
    String fs;
    fs.Printf("curvature.res");
    oss_ss.open(fs.CStr(),std::ios::out);
    oss_ss << Util::_10_6  << "Kappa" << Util::_8s << "Area" << Util::_8s << "Pc" << Util::_8s << "Sr" << std::endl;
    
    for (size_t i=Nmin;i<=Nmax;i++)
    {
        String fname;
        fname.Printf    ("%s_%04d.h5", fileLBM.CStr(), i);

        if (Util::FileExists(fname))
        {
            hid_t file_id;
            file_id = H5Fopen(fname.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
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
            H5Fclose(file_id);
        }
        else throw new Fatal("File <%s> not found",fname.CStr());

        for (size_t j=0;j<Nproc;j++)
        {
            pthread_create(&thrs[j], NULL, LBM::GlobalIni, &MTD[j]);
        }
        for (size_t j=0;j<Nproc;j++)
        {
            pthread_join(thrs[j], NULL);
        }
        for (size_t j=0;j<Nproc;j++)
        {
            pthread_create(&thrs[j], NULL, LBM::GlobalApplyForce, &MTD[j]);
        }
        for (size_t j=0;j<Nproc;j++)
        {
            pthread_join(thrs[j], NULL);
        }

        double kappa = 0.0;
        size_t nc    = 0;

        double Sr = 0.0;

        for (size_t j=0;j<Dom.Lat[0].Cells.Size();j++)
        {
            double wr = Dom.Lat[1].Cells[j]->Rho;
            if (wr>=dthmin)
            {
                Sr+=1.0;
            }
        }

        Sr/=(Dom.Lat[0].Cells.Size()*(1-Dom.Lat[0].SolidFraction()));

        double p1 = 0.0;
        double p2 = 0.0;

        for (size_t i=0;i<Dom.Lat[0].Ndim(1);i++)
        for (size_t j=0;j<Dom.Lat[0].Ndim(2);j++)
        {
            Cell * c0 = Dom.Lat[0].GetCell(iVec3_t(0,i,j));
            Cell * c1 = Dom.Lat[1].GetCell(iVec3_t(0,i,j));

            p1 += (c0->Rho + c1->Rho + Dom.Gmix*(c0->Rho*c1->Rho))/3.0;

            c0 = Dom.Lat[0].GetCell(iVec3_t(Dom.Lat[0].Ndim(0)-1,i,j));
            c1 = Dom.Lat[1].GetCell(iVec3_t(Dom.Lat[0].Ndim(0)-1,i,j));

            p2 += (c0->Rho + c1->Rho + Dom.Gmix*(c0->Rho*c1->Rho))/3.0;
        }

        double Pc = (p1 - p2)/(Dom.Lat[0].Ndim(1)*Dom.Lat[0].Ndim(2));

        double Area = 0.0;

        for (size_t j=0;j<NoSolid.Size();j++)        
        {
            size_t idx = NoSolid[j];
            Cell * c = Dom.Lat[0].Cells[idx];
            if (Dom.Lat[0].Cells[idx]->Rho>=dthmin&&Dom.Lat[0].Cells[idx]->Rho<=dthmax)
            {
                Vec3_t force  = Dom.Lat[0].Cells[idx]->BForce;
                Vec3_t normal = force/norm(force);
                iVec3_t dpos  = ds*normal;
                iVec3_t pos1  = Dom.Lat[0].Cells[idx]->Index + dpos;
                iVec3_t pos2  = Dom.Lat[0].Cells[idx]->Index - dpos;
                if (Dom.Lat[0].GetCell(pos1)->IsSolid) continue;
                if (Dom.Lat[0].GetCell(pos2)->IsSolid) continue;
                double pxmin0 = Dom.Lat[0].GetCell(pos1)->Rho;
                double pxmin1 = Dom.Lat[1].GetCell(pos1)->Rho;
                double pxmax0 = Dom.Lat[0].GetCell(pos2)->Rho;
                double pxmax1 = Dom.Lat[1].GetCell(pos2)->Rho;
                double Pc     = (pxmin0 + pxmin1 - pxmax0 - pxmax1 + Dom.Gmix*(pxmin0*pxmin1 - pxmax0*pxmax1))/3.0;
                kappa += 0.5*Pc/st;
                nc++;
                for (size_t k=1;k<=6;k++)
                {
                    Cell * nb = Dom.Lat[0].Cells[c->Neighs[k]];
                    if (nb->Rho<=dthmin) Area+=1.0;
                }
            }
        }
        //if (nc==0) kappa = 0.0;
        //else       kappa/= nc;

        oss_ss << Util::_8s << kappa << Util::_8s << Area/NoSolid.Size() << Util::_8s << Pc << Util::_8s << Sr << std::endl;
    }


    oss_ss.close();

    return 0;
}
MECHSYS_CATCH


