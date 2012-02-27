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

enum GradCase
{
    Gx,  ///< Gradient in the x direction solely
    Gy,  ///< Gradient in the x direction solely
    Gz   ///< Gradient in the x direction solely
};

inline void VelInter(LBM::Domain & Dom, Vec3_t const & X, Vec3_t & V)
{
    Vec3_t  frac;
    iVec3_t In;
    double dummy;
    int N = 100;
    In(0) = (size_t) X(0);
    In(1) = (size_t) X(1);
    In(2) = (size_t) X(2);
    if ((In(0)>N-2)||(In(1)>N-2)||(In(2)>N-2))
    {
        V = OrthoSys::O;
        return;
    }
    frac(0) = X(0) - In(0);
    frac(1) = X(1) - In(1);
    frac(2) = X(2) - In(2);
    Vec3_t v1,v2,v3,v4,v5,v6;
    v1 = Dom.Lat[0].GetCell(In+iVec3_t(1,0,0))->Vel*frac(0) + Dom.Lat[0].GetCell(In+iVec3_t(0,0,0))->Vel*(1.0-frac(0));
    v2 = Dom.Lat[0].GetCell(In+iVec3_t(1,1,0))->Vel*frac(0) + Dom.Lat[0].GetCell(In+iVec3_t(0,1,0))->Vel*(1.0-frac(0));
    v3 = Dom.Lat[0].GetCell(In+iVec3_t(1,1,1))->Vel*frac(0) + Dom.Lat[0].GetCell(In+iVec3_t(0,1,1))->Vel*(1.0-frac(0));
    v4 = Dom.Lat[0].GetCell(In+iVec3_t(1,0,1))->Vel*frac(0) + Dom.Lat[0].GetCell(In+iVec3_t(0,0,1))->Vel*(1.0-frac(0));
    v5 = v2*frac(1) + v1*(1.0-frac(1));
    v6 = v3*frac(1) + v4*(1.0-frac(1));
    V  = v6*frac(2) + v5*(1.0-frac(2));
}

int main(int argc, char **argv) try
{
    String filename ("graflow_0099.h5");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    int    Case = atoi(argv[1]);
    float  Ndiv = atof(argv[2]);
    float    dt = atof(argv[3]);

    size_t N      = 100;

    LBM::Domain Dom(D3Q15, 0.1, iVec3_t(N,N,N), 1.0, 1.0);
    hid_t file_id;
    file_id = H5Fopen(filename.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);
    float * Density   = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
    float * Gamma     = new float[  Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
    float * Vvec      = new float[3*Dom.Lat[0].Ndim[0]*Dom.Lat[0].Ndim[1]*Dom.Lat[0].Ndim[2]];
    H5LTread_dataset_float(file_id,"Density_0" ,Density);
    H5LTread_dataset_float(file_id,"Gamma_0"   ,Gamma  );
    H5LTread_dataset_float(file_id,"Velocity_0",Vvec   );
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Cell * c = Dom.Lat[0].Cells[i];
        if (Gamma[i]>0.5) c->IsSolid = true;
        c->Gamma = Gamma[i];
        Vec3_t V;
        V(0) =  Vvec[3*i  ];
        V(1) =  Vvec[3*i+1];
        V(2) =  Vvec[3*i+2];
        c->Initialize(Density[i],V);
    }
    Dom.WriteXDMF("test");
    delete [] Density ;
    delete [] Gamma   ;
    delete [] Vvec    ;

    Array <double> Len;
    Array <bool>   Spam;

    if (Case==0)
    {
        for (size_t i=0;i<Ndiv;i++)
        for (size_t j=0;j<Ndiv;j++)
        {
            bool valid = true;
            bool spam  = false;
            double len = 0.0;
            size_t NC  = 0;
            iVec3_t In;
            Vec3_t x(OrthoSys::O);
            x(0) = (0.5 + 0)*(N/Ndiv);
            x(1) = (0.5 + i)*(N/Ndiv);
            x(2) = (0.5 + j)*(N/Ndiv);
            In(0) = (size_t) x(0);
            In(1) = (size_t) x(1);
            In(2) = (size_t) x(2);
            //std::cout << x << In << std::endl;
            Cell * nc = Dom.Lat[0].GetCell(In);
            if (nc->IsSolid) continue;
            while (valid&&(NC<1.0e7))
            {
                Vec3_t  frac;
                double dummy;
                In(0) = (size_t) x(0);
                In(1) = (size_t) x(1);
                In(2) = (size_t) x(2);
                if ((In(1)==0)||(In(1)>=N-2)||(In(2)==0)||(In(2)>=N-2)||(x(0)<0.0)||(x(1)<0.0)||(x(2)<0.0))
                {
                    valid = false;
                    break;
                }
                if (In(0)>=N-1)
                {
                    valid = false;
                    spam  = true;
                    break;
                }
                nc    = Dom.Lat[0].GetCell(In);
                //std::cout << x << In << nc->Vel << " " << len << std::endl;
                if (nc->IsSolid) break;
                //Runge Kutta integration
                Vec3_t k1,k2,k3,k4;
                VelInter(Dom,x       ,k1);
                k1 *= dt;
                VelInter(Dom,x+0.5*k1,k2);
                k2 *= dt;
                VelInter(Dom,x+0.5*k2,k3);
                k3 *= dt;
                VelInter(Dom,x    +k3,k4);
                k4 *= dt;
                x  += (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
                len+= norm((1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4));
                NC++;
            }
            //std::cout << x << " " << len << " " << spam << " " << NC << std::endl;
            Len.Push(len);
            Spam.Push(spam);
        }
    }

    if (Case==1)
    {
        for (size_t i=0;i<Ndiv;i++)
        for (size_t j=0;j<Ndiv;j++)
        {
            bool valid = true;
            bool spam  = false;
            double len = 0.0;
            size_t NC  = 0;
            iVec3_t In;
            Vec3_t x(OrthoSys::O);
            x(0) = (0.5 + i)*(N/Ndiv);
            x(1) = (0.5 + 0)*(N/Ndiv);
            x(2) = (0.5 + j)*(N/Ndiv);
            In(0) = (size_t) x(0);
            In(1) = (size_t) x(1);
            In(2) = (size_t) x(2);
            //std::cout << x << In << std::endl;
            Cell * nc = Dom.Lat[0].GetCell(In);
            if (nc->IsSolid) continue;
            while (valid&&(NC<1.0e7))
            {
                Vec3_t  frac;
                double dummy;
                In(0) = (size_t) x(0);
                In(1) = (size_t) x(1);
                In(2) = (size_t) x(2);
                if ((In(0)==0)||(In(0)>=N-2)||(In(2)==0)||(In(2)>=N-2)||(x(0)<0.0)||(x(1)<0.0)||(x(2)<0.0))
                {
                    valid = false;
                    break;
                }
                if (In(1)>=N-2)
                {
                    valid = false;
                    spam  = true;
                    break;
                }
                nc    = Dom.Lat[0].GetCell(In);
                //std::cout << x << In << nc->Vel << " " << len << std::endl;
                if (nc->IsSolid) break;
                //Runge Kutta integration
                Vec3_t k1,k2,k3,k4;
                VelInter(Dom,x       ,k1);
                k1 *= dt;
                VelInter(Dom,x+0.5*k1,k2);
                k2 *= dt;
                VelInter(Dom,x+0.5*k2,k3);
                k3 *= dt;
                VelInter(Dom,x    +k3,k4);
                k4 *= dt;
                x  += (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
                len+= norm((1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4));
                NC++;
            }
            //std::cout << x << " " << len << " " << spam << " " << NC << std::endl;
            Len.Push(len);
            Spam.Push(spam);
        }
    }

    if (Case==2)
    {
        for (size_t i=0;i<Ndiv;i++)
        for (size_t j=0;j<Ndiv;j++)
        {
            bool valid = true;
            bool spam  = false;
            double len = 0.0;
            size_t NC  = 0;
            iVec3_t In;
            Vec3_t x(OrthoSys::O);
            x(0) = (0.5 + i)*(N/Ndiv);
            x(1) = (0.5 + j)*(N/Ndiv);
            x(2) = (0.5 + 0)*(N/Ndiv);
            In(0) = (size_t) x(0);
            In(1) = (size_t) x(1);
            In(2) = (size_t) x(2);
            //std::cout << x << In << std::endl;
            Cell * nc = Dom.Lat[0].GetCell(In);
            if (nc->IsSolid) continue;
            while (valid&&(NC<1.0e6))
            {
                Vec3_t  frac;
                double dummy;
                In(0) = (size_t) x(0);
                In(1) = (size_t) x(1);
                In(2) = (size_t) x(2);
                if ((In(0)==0)||(In(0)>=N-2)||(In(1)==0)||(In(1)>=N-2)||(x(0)<1.0)||(x(1)<1.0)||(x(2)<0.0))
                {
                    valid = false;
                    break;
                }
                if (In(2)>=N-2)
                {
                    valid = false;
                    spam  = true;
                    break;
                }
                nc    = Dom.Lat[0].GetCell(In);
                //std::cout << x << In << nc->Vel << " " << len << std::endl;
                if (nc->IsSolid) break;
                //Runge Kutta integration
                Vec3_t k1,k2,k3,k4;
                VelInter(Dom,x       ,k1);
                k1 *= dt;
                VelInter(Dom,x+0.5*k1,k2);
                k2 *= dt;
                VelInter(Dom,x+0.5*k2,k3);
                k3 *= dt;
                VelInter(Dom,x    +k3,k4);
                k4 *= dt;
                x  += (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
                len+= norm((1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4));
                NC++;
            }
            std::cout << x << " " << len << " " << spam << " " << NC << std::endl;
            Len.Push(len);
            Spam.Push(spam);
        }
    }
    

    std::ofstream torfile("tortuosity.res");
    torfile << Util::_8s << "Len" << Util::_8s << "Spams?" << std::endl;
    for (size_t i=0;i<Len.Size();i++)
    {
        torfile << Util::_8s << Len[i] << Util::_8s << Spam[i] << std::endl;
    }
    torfile.close();


    return 0;
}
MECHSYS_CATCH


