/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;


struct UserData
{
    std::ofstream oss_st;       ///< file for surface tension calculation
    double           rho;
};

void Report (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double max1 = 0.0;
    double min1 = (int) dom.Lat[0].Ndim(0);
    iVec3_t Nd  = dom.Lat[0].Ndim;

    for (size_t i=0;i<Nd(0)/2;i++)
    {
        Cell * c1 = dom.Lat[0].GetCell(iVec3_t(i,Nd(1)/2,0));;
        if (c1->Rho>0.5*dat.rho) 
        {
            max1 = i;
            break;
        }
    }
    for (size_t i=0;i<Nd(0)/2;i++)
    {
        Cell * c1 = dom.Lat[0].GetCell(iVec3_t(Nd(0)-i-1,Nd(1)/2,0));;
        if (c1->Rho>0.5*dat.rho) 
        {
            min1 = Nd(0)-i-1;
            break;
        }
    }

    double Pin = (dom.Lat[0].GetCell(Nd/2          )->Rho + dom.Lat[1].GetCell(Nd/2          )->Rho
                 + dom.Gmix*dom.Lat[0].GetCell(Nd/2          )->Rho*dom.Lat[1].GetCell(Nd/2)->Rho          )/3.0;
    double Pout= (dom.Lat[0].GetCell(iVec3_t(0,0,0))->Rho + dom.Lat[1].GetCell(iVec3_t(0,0,0))->Rho
                 + dom.Gmix*dom.Lat[0].GetCell(iVec3_t(0,0,0))->Rho*dom.Lat[1].GetCell(iVec3_t(0,0,0))->Rho)/3.0;
    double R   = 0.5*(min1 - max1);
    
    dat.oss_st << dom.Time << Util::_8s << Pin << Util::_8s << Pout << Util::_8s << R << Util::_8s << R*(Pin-Pout) << std::endl;
}



int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    double Gmix  = 0.001;
    double R     = 10.0;
    double rho   = 1.0;
    double rhodif= 0.01;
    double Tf    = 5000.0;
    double dtout = 50.0;
    if (argc>=2)
    {
        Nproc  =atoi(argv[1]);
        Gmix   =atof(argv[2]);
        R      =atof(argv[3]);
        rho    =atof(argv[4]);
        rhodif =atof(argv[5]);
        Tf     =atof(argv[6]);
        dtout  =atof(argv[7]);
    }
    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    size_t nx = 128, ny = 128;

    // Setting top and bottom wall as solid
    LBM::Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.rho      = rho;
    int obsX = nx/2, obsY = ny/2;
    int radius =  (int) R;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(rho-rhodif,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(rhodif    ,V);
		}
		else
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(rhodif    ,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(rho-rhodif,V);
		}
    }

    // Set parameters
    Dom.Gmix     = Gmix;
    Dom.Lat[0].G = 0.0;
    Dom.Lat[1].G = 0.0;

    String fs;
    fs.Printf("surface_tension.res");
    dat.oss_st.open(fs.CStr(),std::ios::out);
    dat.oss_st << Util::_10_6  <<  "Time" << Util::_8s << "Pin" << Util::_8s << "Pout" << Util::_8s << "R" << Util::_8s << "ST" << std::endl;
    Dom.Solve(Tf,dtout,NULL,Report,"multicomp",true,Nproc);
    dat.oss_st.close();

    return 0;
}
MECHSYS_CATCH
